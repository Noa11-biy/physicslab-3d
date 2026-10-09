#include "physicslab/render/Renderer.hpp"

#include <glad/gl.h>

#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <sstream>

namespace pl {
namespace {

bool readFile(const std::string& path, std::string& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "Renderer : impossible de lire %s\n", path.c_str());
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    out = ss.str();
    return true;
}

GLuint compileShader(GLenum type, const std::string& path) {
    std::string src;
    if (!readFile(path, src)) return 0;

    const GLuint shader = glCreateShader(type);
    const char* text = src.c_str();
    glShaderSource(shader, 1, &text, nullptr);
    glCompileShader(shader);

    GLint ok = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        std::fprintf(stderr, "Renderer : erreur de compilation de %s :\n%s\n", path.c_str(), log);
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

}  // namespace

bool Renderer::init(const std::string& shaderDir) {
    const GLuint vs = compileShader(GL_VERTEX_SHADER, shaderDir + "/line.vert");
    const GLuint fs = compileShader(GL_FRAGMENT_SHADER, shaderDir + "/line.frag");
    if (!vs || !fs) return false;

    program_ = glCreateProgram();
    glAttachShader(program_, vs);
    glAttachShader(program_, fs);
    glLinkProgram(program_);
    glDeleteShader(vs);
    glDeleteShader(fs);

    GLint ok = GL_FALSE;
    glGetProgramiv(program_, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetProgramInfoLog(program_, sizeof(log), nullptr, log);
        std::fprintf(stderr, "Renderer : erreur d'édition de liens :\n%s\n", log);
        return false;
    }

    uViewProj_ = glGetUniformLocation(program_, "uViewProj");
    uPointSize_ = glGetUniformLocation(program_, "uPointSize");
    uRound_ = glGetUniformLocation(program_, "uRound");

    // État direct (DSA, OpenGL 4.5) : un VAO et un VBO dynamique partagés par tous les dessins.
    glCreateVertexArrays(1, &vao_);
    glCreateBuffers(1, &vbo_);
    glVertexArrayVertexBuffer(vao_, 0, vbo_, 0, sizeof(Vertex));

    glEnableVertexArrayAttrib(vao_, 0);
    glVertexArrayAttribFormat(vao_, 0, 3, GL_FLOAT, GL_FALSE, static_cast<GLuint>(offsetof(Vertex, x)));
    glVertexArrayAttribBinding(vao_, 0, 0);

    glEnableVertexArrayAttrib(vao_, 1);
    glVertexArrayAttribFormat(vao_, 1, 3, GL_FLOAT, GL_FALSE, static_cast<GLuint>(offsetof(Vertex, r)));
    glVertexArrayAttribBinding(vao_, 1, 0);
    return true;
}

void Renderer::shutdown() {
    if (vbo_) glDeleteBuffers(1, &vbo_);
    if (vao_) glDeleteVertexArrays(1, &vao_);
    if (program_) glDeleteProgram(program_);
    vbo_ = vao_ = program_ = 0;
    vboCapacity_ = vboOffset_ = 0;
}

void Renderer::beginFrame(int fbW, int fbH, const ViewRect& view, const Camera& camera,
                          const float clearColor[3]) {
    // 1) Efface toute la fenêtre.
    glDisable(GL_SCISSOR_TEST);
    glViewport(0, 0, fbW, fbH);
    glClearColor(clearColor[0], clearColor[1], clearColor[2], 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // 2) Restreint le dessin 3D à la zone libre (hors panneaux). Origine OpenGL en bas à gauche.
    const int glY = fbH - view.y - view.h;
    glViewport(view.x, glY, view.w, view.h);
    glEnable(GL_SCISSOR_TEST);
    glScissor(view.x, glY, view.w, view.h);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);  // à profondeur égale, le dernier dessin gagne (courbes superposées)
    glEnable(GL_PROGRAM_POINT_SIZE);

    float viewProj[16];
    const float aspect = view.h > 0 ? static_cast<float>(view.w) / static_cast<float>(view.h) : 1.0f;
    camera.viewProjection(aspect, viewProj);

    glUseProgram(program_);
    glUniformMatrix4fv(uViewProj_, 1, GL_FALSE, viewProj);
    glBindVertexArray(vao_);
    vboOffset_ = 0;
}

void Renderer::draw(Primitive primitive, const std::vector<Vertex>& vertices, float pointSize) {
    if (vertices.empty()) return;

    const long long bytes = static_cast<long long>(vertices.size() * sizeof(Vertex));
    if (vboOffset_ + bytes > vboCapacity_) {  // plus de place dans cette image : nouvelle allocation (4 Mio au moins)
        vboCapacity_ = std::max<long long>({4LL << 20, 2 * vboCapacity_, bytes});
        glNamedBufferData(vbo_, static_cast<GLsizeiptr>(vboCapacity_), nullptr, GL_STREAM_DRAW);
        vboOffset_ = 0;
    }
    glNamedBufferSubData(vbo_, static_cast<GLintptr>(vboOffset_), static_cast<GLsizeiptr>(bytes), vertices.data());
    glVertexArrayVertexBuffer(vao_, 0, vbo_, static_cast<GLintptr>(vboOffset_), sizeof(Vertex));
    vboOffset_ += bytes;
    glUniform1f(uPointSize_, pointSize);
    glUniform1i(uRound_, primitive == Primitive::Points ? 1 : 0);

    GLenum mode = GL_LINES;
    if (primitive == Primitive::LineStrip) mode = GL_LINE_STRIP;
    if (primitive == Primitive::Points) mode = GL_POINTS;
    glDrawArrays(mode, 0, static_cast<GLsizei>(vertices.size()));
}

std::vector<Vertex> makeGrid(int halfCells, float cellSize) {
    std::vector<Vertex> v;
    const float extent = static_cast<float>(halfCells) * cellSize;
    const float shade = 0.22f;
    for (int i = -halfCells; i <= halfCells; ++i) {
        const float p = static_cast<float>(i) * cellSize;
        v.push_back({p, 0.0f, -extent, shade, shade, shade});  // ligne parallèle à Z
        v.push_back({p, 0.0f, extent, shade, shade, shade});
        v.push_back({-extent, 0.0f, p, shade, shade, shade});  // ligne parallèle à X
        v.push_back({extent, 0.0f, p, shade, shade, shade});
    }
    return v;
}

std::vector<Vertex> makeAxes(float length) {
    return {{0, 0, 0, 0.90f, 0.25f, 0.25f}, {length, 0, 0, 0.90f, 0.25f, 0.25f},
            {0, 0, 0, 0.30f, 0.80f, 0.30f}, {0, length, 0, 0.30f, 0.80f, 0.30f},
            {0, 0, 0, 0.30f, 0.45f, 0.95f}, {0, 0, length, 0.30f, 0.45f, 0.95f}};
}

}  // namespace pl
