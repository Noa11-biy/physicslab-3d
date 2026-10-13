#include "physicslab/render/GpuWave2D.hpp"

#include <glad/gl.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "physicslab/render/GpuProgram.hpp"

namespace pl {
namespace {

constexpr int kMaxSources = 4;  // taille des tableaux uniformes du shader
constexpr int kMaxLines = 2;

std::vector<float> toFloat(const std::vector<double>& v) { return std::vector<float>(v.begin(), v.end()); }

void fill(unsigned int& buffer, const std::vector<float>& data) {
    if (buffer == 0) glCreateBuffers(1, &buffer);
    glNamedBufferData(buffer, static_cast<GLsizeiptr>(data.size() * sizeof(float)), data.data(), GL_DYNAMIC_DRAW);
}

}  // namespace

bool GpuWave2D::init(const std::string& shaderDir) {
    shutdown();
    program_ = buildComputeProgram(shaderDir + "/wave2d.comp", "", "GpuWave2D");
    if (program_ == 0) return false;
    uNx_ = glGetUniformLocation(program_, "uNx");
    uNy_ = glGetUniformLocation(program_, "uNy");
    uDt_ = glGetUniformLocation(program_, "uDt");
    uFree_ = glGetUniformLocation(program_, "uFree");
    uSourceCount_ = glGetUniformLocation(program_, "uSourceCount");
    uSourceNode_ = glGetUniformLocation(program_, "uSourceNode[0]");
    uSourceValue_ = glGetUniformLocation(program_, "uSourceValue[0]");
    uLineCount_ = glGetUniformLocation(program_, "uLineCount");
    uLineNode_ = glGetUniformLocation(program_, "uLineNode[0]");
    uLineValue_ = glGetUniformLocation(program_, "uLineValue[0]");
    return true;
}

void GpuWave2D::shutdown() {
    for (unsigned int& b : field_) {
        if (b) glDeleteBuffers(1, &b);
        b = 0;
    }
    for (unsigned int* b : {&coefX_, &coefY_, &sigma_}) {
        if (*b) glDeleteBuffers(1, b);
        *b = 0;
    }
    if (program_) glDeleteProgram(program_);
    program_ = 0;
}

void GpuWave2D::upload(const waves::Wave2D& wave) {
    nx_ = wave.pointsX();
    ny_ = wave.pointsY();
    cellsX_ = nx_ - 1;
    dx_ = wave.dx();
    dt_ = wave.dt();
    ratio_ = dt_ * dt_ / (dx_ * dx_);
    steps_ = wave.steps();
    const waves::Wave2DParams& p = wave.params();
    freeEdge_[0] = p.left == waves::Edge::Free;
    freeEdge_[1] = p.right == waves::Edge::Free;
    freeEdge_[2] = p.bottom == waves::Edge::Free;
    freeEdge_[3] = p.top == waves::Edge::Free;
    sources_ = wave.sources();
    lines_ = wave.lineSources();
    params_ = p;
    params_.speedField.clear();  // inutile ici (les coefficients de faces suffisent) : pas de copie de plus
    params_.wall.clear();

    fill(field_[0], toFloat(wave.uPrevious()));
    fill(field_[1], toFloat(wave.u()));
    fill(field_[2], std::vector<float>(wave.u().size(), 0.0f));
    fill(coefX_, toFloat(wave.coefficientX()));
    fill(coefY_, toFloat(wave.coefficientY()));
    fill(sigma_, toFloat(wave.sigmaField()));
    slot_ = 1;
}

// Valeur injectée à chaque source à l'instant courant : ratio * A * rampe(t) * sin(w t), calculée en double comme Wave2D::step ; noeud le plus proche,
// jamais sur un bord fixe.
void GpuWave2D::sourceValues(int* nodes, float* values, int& count) const {
    count = 0;
    const double t = steps_ * dt_;
    const double pi = 3.14159265358979323846;
    for (const waves::PointSource& src : sources_) {
        if (count >= kMaxSources) break;
        const int i = std::clamp(static_cast<int>(std::lround(src.x / dx_)), 0, nx_ - 1), j = std::clamp(static_cast<int>(std::lround(src.y / dx_)), 0, ny_ - 1);
        const bool fixedNode = (i == 0 && !freeEdge_[0]) || (i == nx_ - 1 && !freeEdge_[1]) || (j == 0 && !freeEdge_[2]) || (j == ny_ - 1 && !freeEdge_[3]);
        if (fixedNode) continue;
        const double period = 2.0 * pi / src.omega;
        const double r = src.rampPeriods > 0.0 ? std::min(1.0, t / (src.rampPeriods * period)) : 1.0;
        const double ramp = std::sin(0.5 * pi * r) * std::sin(0.5 * pi * r);
        nodes[2 * count] = i;
        nodes[2 * count + 1] = j;
        values[count] = static_cast<float>(ratio_ * src.amplitude * ramp * std::sin(src.omega * t));
        ++count;
    }
}

// Sources en ligne (O4) : (colonne, j0, j1) et valeur (dt² / dx) A rampe(t) sin(w t) = ratio * dx * ..., en double comme Wave2D::step.
void GpuWave2D::lineValues(int* nodes, float* values, int& count) const {
    count = 0;
    const double t = steps_ * dt_;
    for (const waves::LineSource& src : lines_) {
        if (count >= kMaxLines) break;
        int i, j0, j1;
        if (!waves::lineSourceNodes(src, params_, i, j0, j1)) continue;
        nodes[3 * count] = i;
        nodes[3 * count + 1] = j0;
        nodes[3 * count + 2] = j1;
        values[count] = static_cast<float>(ratio_ * dx_ * src.amplitude * waves::sourceRamp(t, src.omega, src.rampPeriods) * std::sin(src.omega * t));
        ++count;
    }
}

void GpuWave2D::bindBuffers() const {
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, field_[(slot_ + 2) % 3]);  // précédent
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, field_[slot_]);            // courant
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, field_[(slot_ + 1) % 3]);  // suivant
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, coefX_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, coefY_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, sigma_);
}

void GpuWave2D::step(int count) {
    if (!ready() || nx_ == 0) return;
    glUseProgram(program_);
    glUniform1i(uNx_, nx_);
    glUniform1i(uNy_, ny_);
    glUniform1f(uDt_, static_cast<float>(dt_));
    glUniform4i(uFree_, freeEdge_[0], freeEdge_[1], freeEdge_[2], freeEdge_[3]);
    const GLuint groupsX = static_cast<GLuint>((nx_ + 15) / 16), groupsY = static_cast<GLuint>((ny_ + 15) / 16);
    for (int n = 0; n < count; ++n) {
        int nodes[2 * kMaxSources] = {0};
        float values[kMaxSources] = {0.0f};
        int sourceCount = 0;
        sourceValues(nodes, values, sourceCount);
        glUniform1i(uSourceCount_, sourceCount);
        if (sourceCount > 0) {
            glUniform2iv(uSourceNode_, sourceCount, nodes);
            glUniform1fv(uSourceValue_, sourceCount, values);
        }
        int lineNodes[3 * kMaxLines] = {0};
        float lineVals[kMaxLines] = {0.0f};
        int lineCount = 0;
        lineValues(lineNodes, lineVals, lineCount);
        glUniform1i(uLineCount_, lineCount);
        if (lineCount > 0) {
            glUniform3iv(uLineNode_, lineCount, lineNodes);
            glUniform1fv(uLineValue_, lineCount, lineVals);
        }
        bindBuffers();
        glDispatchCompute(groupsX, groupsY, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);  // le champ écrit ici est lu par le pas suivant, puis relu par l'hôte
        slot_ = (slot_ + 1) % 3;
        ++steps_;
    }
}

void GpuWave2D::readField(std::vector<float>& out) const {
    out.resize(static_cast<std::size_t>(nx_) * static_cast<std::size_t>(ny_));
    glGetNamedBufferSubData(field_[slot_], 0, static_cast<GLsizeiptr>(out.size() * sizeof(float)), out.data());
}

double GpuWave2D::benchmark(int steps) {
    GLuint query = 0;
    glGenQueries(1, &query);
    glBeginQuery(GL_TIME_ELAPSED, query);
    step(steps);
    glEndQuery(GL_TIME_ELAPSED);
    GLuint64 nanoseconds = 0;
    glGetQueryObjectui64v(query, GL_QUERY_RESULT, &nanoseconds);  // attend la fin
    glDeleteQueries(1, &query);
    return static_cast<double>(nanoseconds) * 1e-9;
}

}  // namespace pl
