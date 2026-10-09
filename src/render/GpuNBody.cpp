#include "physicslab/render/GpuNBody.hpp"

#include <glad/gl.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

namespace pl {
namespace {

// Source du compute shader, complétée juste après la ligne #version par les #define de configuration.
bool loadShaderSource(const std::string& path, const std::string& defines, std::string& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "GpuNBody : impossible de lire %s\n", path.c_str());
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    const std::size_t eol = text.find('\n');
    if (text.compare(0, 8, "#version") != 0 || eol == std::string::npos) {
        std::fprintf(stderr, "GpuNBody : %s doit commencer par #version\n", path.c_str());
        return false;
    }
    out = text.substr(0, eol + 1) + defines + text.substr(eol + 1);
    return true;
}

// Compile et lie un programme de calcul. Renvoie 0 en cas d'échec (message sur stderr).
GLuint buildComputeProgram(const std::string& path, const std::string& defines) {
    std::string source;
    if (!loadShaderSource(path, defines, source)) return 0;

    const GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    const char* text = source.c_str();
    glShaderSource(shader, 1, &text, nullptr);
    glCompileShader(shader);
    GLint ok = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        std::fprintf(stderr, "GpuNBody : erreur de compilation de %s :\n%s\n", path.c_str(), log);
        glDeleteShader(shader);
        return 0;
    }

    const GLuint program = glCreateProgram();
    glAttachShader(program, shader);
    glLinkProgram(program);
    glDeleteShader(shader);
    glGetProgramiv(program, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetProgramInfoLog(program, sizeof(log), nullptr, log);
        std::fprintf(stderr, "GpuNBody : erreur d'édition de liens de %s :\n%s\n", path.c_str(), log);
        glDeleteProgram(program);
        return 0;
    }
    return program;
}

}  // namespace

GpuLimits queryGpuLimits() {
    GpuLimits g;
    const auto text = [](GLenum name) {
        const GLubyte* s = glGetString(name);
        return std::string(s ? reinterpret_cast<const char*>(s) : "?");
    };
    g.renderer = text(GL_RENDERER);
    g.version = text(GL_VERSION);

    GLint v = 0;
    glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &v);
    g.maxGroupInvocations = v;
    for (int axis = 0; axis < 3; ++axis) {
        v = 0;
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, static_cast<GLuint>(axis), &v);
        g.maxGroupSize[axis] = v;
    }
    v = 0;
    glGetIntegerv(GL_MAX_COMPUTE_SHARED_MEMORY_SIZE, &v);
    g.maxSharedBytes = v;
    GLint64 big = 0;
    glGetInteger64v(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &big);
    g.maxStorageBlockBytes = big;

    GLint count = 0;
    glGetIntegerv(GL_NUM_EXTENSIONS, &count);
    for (GLint i = 0; i < count; ++i) {
        const char* ext = reinterpret_cast<const char*>(glGetStringi(GL_EXTENSIONS, static_cast<GLuint>(i)));
        if (ext && std::strcmp(ext, "GL_ARB_gpu_shader_fp64") == 0) g.doublePrecision = true;
    }
    return g;
}

void centerOnBarycenter(const double* positions, const double* masses, int n, double* out, double* center) {
    double cx = 0.0, cy = 0.0, cz = 0.0, total = 0.0;
    for (int i = 0; i < n; ++i) {
        cx += masses[i] * positions[3 * i];
        cy += masses[i] * positions[3 * i + 1];
        cz += masses[i] * positions[3 * i + 2];
        total += masses[i];
    }
    if (total > 0.0) {
        cx /= total; cy /= total; cz /= total;
    } else if (n > 0) {  // masses nulles : moyenne simple
        for (int i = 0; i < n; ++i) { cx += positions[3 * i]; cy += positions[3 * i + 1]; cz += positions[3 * i + 2]; }
        cx /= n; cy /= n; cz /= n;
    }
    for (int i = 0; i < n; ++i) {
        out[3 * i] = positions[3 * i] - cx;
        out[3 * i + 1] = positions[3 * i + 1] - cy;
        out[3 * i + 2] = positions[3 * i + 2] - cz;
    }
    if (center) {
        center[0] = cx;
        center[1] = cy;
        center[2] = cz;
    }
}

bool GpuNBody::init(const std::string& shaderDir, GpuPrecision precision, int workGroupSize) {
    shutdown();
    precision_ = precision;
    wg_ = workGroupSize;

    std::string defines = "#define WG " + std::to_string(wg_) + "\n";
    if (precision_ == GpuPrecision::Double) defines += "#define GPU_DOUBLE\n";
    program_ = buildComputeProgram(shaderDir + "/nbody.comp", defines);
    stepProgram_ = buildComputeProgram(shaderDir + "/nbody_step.comp", defines);
    if (!program_ || !stepProgram_) {
        shutdown();
        return false;
    }

    uCount_ = glGetUniformLocation(program_, "uCount");
    uFirst_ = glGetUniformLocation(program_, "uFirstGroup");
    uG_ = glGetUniformLocation(program_, "uG");
    uEps2_ = glGetUniformLocation(program_, "uEps2");
    uStepCount_ = glGetUniformLocation(stepProgram_, "uCount");
    uDtKick_ = glGetUniformLocation(stepProgram_, "uDtKick");
    uDtDrift_ = glGetUniformLocation(stepProgram_, "uDtDrift");
    return true;
}

void GpuNBody::shutdown() {
    if (bodies_) glDeleteBuffers(1, &bodies_);
    if (accel_) glDeleteBuffers(1, &accel_);
    if (velocities_) glDeleteBuffers(1, &velocities_);
    if (!queryPool_.empty()) glDeleteQueries(static_cast<GLsizei>(queryPool_.size()), queryPool_.data());
    if (program_) glDeleteProgram(program_);
    if (stepProgram_) glDeleteProgram(stepProgram_);
    bodies_ = accel_ = velocities_ = program_ = stepProgram_ = 0;
    queryPool_.clear();
    pending_.clear();
    queuedSeconds_ = 0.0;
    n_ = capacity_ = 0;
    stateReady_ = accelValid_ = false;
}

void GpuNBody::setReal(unsigned int program, int location, double value) const {
    if (precision_ == GpuPrecision::Double) {
        glProgramUniform1d(program, location, value);
    } else {
        glProgramUniform1f(program, location, static_cast<float>(value));
    }
}

void GpuNBody::uploadBodies(const double* positions, const double* masses, int n) {
    n_ = n;
    accelValid_ = false;
    stateReady_ = false;
    if (!program_ || n <= 0) return;

    const std::size_t bytesPerBody = (precision_ == GpuPrecision::Double ? 4 * sizeof(double) : 4 * sizeof(float));
    if (n > capacity_) {
        if (bodies_) glDeleteBuffers(1, &bodies_);
        if (accel_) glDeleteBuffers(1, &accel_);
        if (velocities_) glDeleteBuffers(1, &velocities_);
        glCreateBuffers(1, &bodies_);
        glCreateBuffers(1, &accel_);
        glCreateBuffers(1, &velocities_);
        const GLsizeiptr bytes = static_cast<GLsizeiptr>(bytesPerBody * static_cast<std::size_t>(n));
        glNamedBufferData(bodies_, bytes, nullptr, GL_DYNAMIC_DRAW);
        glNamedBufferData(accel_, bytes, nullptr, GL_DYNAMIC_DRAW);
        glNamedBufferData(velocities_, bytes, nullptr, GL_DYNAMIC_DRAW);
        capacity_ = n;
    }

    std::vector<double> centered(3 * static_cast<std::size_t>(n));
    offset_[0] = offset_[1] = offset_[2] = 0.0;
    if (center_) {
        centerOnBarycenter(positions, masses, n, centered.data(), offset_);
    } else {
        std::copy(positions, positions + 3 * n, centered.begin());
    }

    if (precision_ == GpuPrecision::Double) {
        std::vector<double> host(4 * static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            host[4 * i] = centered[3 * i];
            host[4 * i + 1] = centered[3 * i + 1];
            host[4 * i + 2] = centered[3 * i + 2];
            host[4 * i + 3] = masses[i];
        }
        glNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(double)), host.data());
    } else {
        std::vector<float> host(4 * static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            host[4 * i] = static_cast<float>(centered[3 * i]);
            host[4 * i + 1] = static_cast<float>(centered[3 * i + 1]);
            host[4 * i + 2] = static_cast<float>(centered[3 * i + 2]);
            host[4 * i + 3] = static_cast<float>(masses[i]);
        }
        glNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(float)), host.data());
    }
}

void GpuNBody::setBodies(const double* positions, const double* masses, int n) { uploadBodies(positions, masses, n); }

void GpuNBody::setState(const double* positions, const double* velocities, const double* masses, int n) {
    uploadBodies(positions, masses, n);
    if (!program_ || n <= 0) return;

    if (precision_ == GpuPrecision::Double) {
        std::vector<double> host(4 * static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
            host[4 * i] = velocities[3 * i];
            host[4 * i + 1] = velocities[3 * i + 1];
            host[4 * i + 2] = velocities[3 * i + 2];
        }
        glNamedBufferSubData(velocities_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(double)), host.data());
    } else {
        std::vector<float> host(4 * static_cast<std::size_t>(n), 0.0f);
        for (int i = 0; i < n; ++i) {
            host[4 * i] = static_cast<float>(velocities[3 * i]);
            host[4 * i + 1] = static_cast<float>(velocities[3 * i + 1]);
            host[4 * i + 2] = static_cast<float>(velocities[3 * i + 2]);
        }
        glNamedBufferSubData(velocities_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(float)), host.data());
    }
    stateReady_ = true;
}

void GpuNBody::dispatchAccel(double G, double softening, bool blocking) {
    glUseProgram(program_);
    glProgramUniform1i(program_, uCount_, n_);
    setReal(program_, uG_, G);
    setReal(program_, uEps2_, softening * softening);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, bodies_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, accel_);

    const int groups = (n_ + wg_ - 1) / wg_;
    const double perGroup = static_cast<double>(wg_) * static_cast<double>(n_);  // interactions d'un groupe
    for (int first = 0; first < groups;) {
        int count = 1;  // débit inconnu : un seul groupe pour le mesurer
        if (rate_ > 0.0) count = std::clamp(static_cast<int>(rate_ * targetDispatchSeconds_ / perGroup), 1, groups - first);

        if (pending_.size() == queryPool_.size()) {
            GLuint q = 0;
            glGenQueries(1, &q);
            queryPool_.push_back(q);
        }
        const GLuint query = queryPool_[pending_.size()];
        glProgramUniform1ui(program_, uFirst_, static_cast<GLuint>(first));
        glBeginQuery(GL_TIME_ELAPSED, query);
        glDispatchCompute(static_cast<GLuint>(count), 1, 1);
        glEndQuery(GL_TIME_ELAPSED);
        pending_.emplace_back(query, count * perGroup);
        queuedSeconds_ += rate_ > 0.0 ? count * perGroup / rate_ : 0.0;
        first += count;

        // Débit inconnu : on attend pour le mesurer avant de dimensionner les envois suivants. Sinon on n'attend que si le travail
        // empilé devient long (pilote réactif, requêtes en nombre raisonnable).
        if (blocking || rate_ <= 0.0 || queuedSeconds_ > 0.25 || pending_.size() >= 256) collectTimings();
    }
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT | GL_BUFFER_UPDATE_BARRIER_BIT);  // accel_ est lu ensuite par d'autres envois
}

void GpuNBody::collectTimings() {
    double interactions = 0.0, seconds = 0.0;
    for (const auto& [query, work] : pending_) {
        GLuint64 nanoseconds = 0;
        glGetQueryObjectui64v(query, GL_QUERY_RESULT, &nanoseconds);  // attend la fin de cet envoi
        seconds += static_cast<double>(nanoseconds) * 1e-9;
        interactions += work;
    }
    kernelSeconds_ += seconds;
    if (seconds > 0.0) rate_ = interactions / seconds;
    pending_.clear();
    queuedSeconds_ = 0.0;
}

void GpuNBody::compute(double G, double softening) {
    kernelSeconds_ = 0.0;
    if (!program_ || n_ <= 0) return;
    dispatchAccel(G, softening, true);
    collectTimings();
    accelValid_ = true;
    accelG_ = G;
    accelEps_ = softening;
}

void GpuNBody::dispatchStepKernel(double kick, double drift) {
    glUseProgram(stepProgram_);
    glProgramUniform1i(stepProgram_, uStepCount_, n_);
    setReal(stepProgram_, uDtKick_, kick);
    setReal(stepProgram_, uDtDrift_, drift);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, bodies_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, accel_);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, velocities_);
    glDispatchCompute(static_cast<GLuint>((n_ + wg_ - 1) / wg_), 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);  // les positions écrites ici sont lues par le calcul des accélérations
}

void GpuNBody::step(double G, double softening, double dt, int steps) {
    kernelSeconds_ = 0.0;
    if (!ready() || !stateReady_ || n_ <= 0 || steps <= 0) return;

    if (!accelValid_ || accelG_ != G || accelEps_ != softening) dispatchAccel(G, softening, true);  // a(x0)
    for (int s = 0; s < steps; ++s) {
        // Deux demi-coups consécutifs (fin du pas précédent, début de celui-ci) fusionnent en un coup entier.
        dispatchStepKernel(s == 0 ? 0.5 * dt : dt, dt);
        dispatchAccel(G, softening, false);  // a(x_{s+1}) ; les envois s'empilent, sans attendre le GPU à chaque pas
    }
    dispatchStepKernel(0.5 * dt, 0.0);  // dernier demi-coup : l'état est synchrone
    collectTimings();                   // attend la fin du lot
    accelValid_ = true;
    accelG_ = G;
    accelEps_ = softening;
}

void GpuNBody::readAccelerations(double* acc) const {
    if (!program_ || n_ <= 0) return;
    if (precision_ == GpuPrecision::Double) {
        std::vector<double> host(4 * static_cast<std::size_t>(n_));
        glGetNamedBufferSubData(accel_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(double)), host.data());
        for (int i = 0; i < n_; ++i) {
            acc[3 * i] = host[4 * i];
            acc[3 * i + 1] = host[4 * i + 1];
            acc[3 * i + 2] = host[4 * i + 2];
        }
    } else {
        std::vector<float> host(4 * static_cast<std::size_t>(n_));
        glGetNamedBufferSubData(accel_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(float)), host.data());
        for (int i = 0; i < n_; ++i) {
            acc[3 * i] = static_cast<double>(host[4 * i]);
            acc[3 * i + 1] = static_cast<double>(host[4 * i + 1]);
            acc[3 * i + 2] = static_cast<double>(host[4 * i + 2]);
        }
    }
}

void GpuNBody::readState(double* positions, double* velocities) const {
    if (!program_ || n_ <= 0) return;
    if (precision_ == GpuPrecision::Double) {
        std::vector<double> host(4 * static_cast<std::size_t>(n_));
        glGetNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(double)), host.data());
        for (int i = 0; i < n_; ++i)
            for (int k = 0; k < 3; ++k) positions[3 * i + k] = host[4 * i + k] + offset_[k];
        glGetNamedBufferSubData(velocities_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(double)), host.data());
        for (int i = 0; i < n_; ++i)
            for (int k = 0; k < 3; ++k) velocities[3 * i + k] = host[4 * i + k];
    } else {
        std::vector<float> host(4 * static_cast<std::size_t>(n_));
        glGetNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(float)), host.data());
        for (int i = 0; i < n_; ++i)
            for (int k = 0; k < 3; ++k) positions[3 * i + k] = static_cast<double>(host[4 * i + k]) + offset_[k];
        glGetNamedBufferSubData(velocities_, 0, static_cast<GLsizeiptr>(host.size() * sizeof(float)), host.data());
        for (int i = 0; i < n_; ++i)
            for (int k = 0; k < 3; ++k) velocities[3 * i + k] = static_cast<double>(host[4 * i + k]);
    }
}

double GpuNBody::potentialEnergy() const {
    if (!program_ || n_ <= 0 || !accelValid_) return 0.0;
    double u = 0.0;
    if (precision_ == GpuPrecision::Double) {
        std::vector<double> acc(4 * static_cast<std::size_t>(n_)), body(acc.size());
        glGetNamedBufferSubData(accel_, 0, static_cast<GLsizeiptr>(acc.size() * sizeof(double)), acc.data());
        glGetNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(body.size() * sizeof(double)), body.data());
        for (int i = 0; i < n_; ++i) u += body[4 * i + 3] * acc[4 * i + 3];
    } else {
        std::vector<float> acc(4 * static_cast<std::size_t>(n_)), body(acc.size());
        glGetNamedBufferSubData(accel_, 0, static_cast<GLsizeiptr>(acc.size() * sizeof(float)), acc.data());
        glGetNamedBufferSubData(bodies_, 0, static_cast<GLsizeiptr>(body.size() * sizeof(float)), body.data());
        for (int i = 0; i < n_; ++i) u += static_cast<double>(body[4 * i + 3]) * static_cast<double>(acc[4 * i + 3]);
    }
    return 0.5 * u;
}

void GpuNBody::accelerations(const double* positions, const double* masses, int n, double G, double softening, double* acc) {
    setBodies(positions, masses, n);
    compute(G, softening);
    readAccelerations(acc);
}

}  // namespace pl
