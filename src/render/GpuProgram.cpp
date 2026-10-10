#include "physicslab/render/GpuProgram.hpp"

#include <glad/gl.h>

#include <cstdio>
#include <fstream>
#include <sstream>

namespace pl {
namespace {

// Source du compute shader, complétée juste après la ligne #version par les #define de configuration.
bool loadShaderSource(const std::string& path, const std::string& defines, const char* tag, std::string& out) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        std::fprintf(stderr, "%s : impossible de lire %s\n", tag, path.c_str());
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    const std::size_t eol = text.find('\n');
    if (text.compare(0, 8, "#version") != 0 || eol == std::string::npos) {
        std::fprintf(stderr, "%s : %s doit commencer par #version\n", tag, path.c_str());
        return false;
    }
    out = text.substr(0, eol + 1) + defines + text.substr(eol + 1);
    return true;
}

}  // namespace

unsigned int buildComputeProgram(const std::string& path, const std::string& defines, const char* tag) {
    std::string source;
    if (!loadShaderSource(path, defines, tag, source)) return 0;

    const GLuint shader = glCreateShader(GL_COMPUTE_SHADER);
    const char* text = source.c_str();
    glShaderSource(shader, 1, &text, nullptr);
    glCompileShader(shader);
    GLint ok = GL_FALSE;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[2048];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        std::fprintf(stderr, "%s : erreur de compilation de %s :\n%s\n", tag, path.c_str(), log);
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
        std::fprintf(stderr, "%s : erreur d'édition de liens de %s :\n%s\n", tag, path.c_str(), log);
        glDeleteProgram(program);
        return 0;
    }
    return program;
}

}  // namespace pl
