#include "physicslab/render/Camera.hpp"

#include <algorithm>
#include <cmath>

namespace pl {
namespace {

struct F3 {
    float x, y, z;
};
F3 operator-(F3 a, F3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
float dot(F3 a, F3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
F3 cross(F3 a, F3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
F3 normalize(F3 a) {
    const float n = std::sqrt(dot(a, a));
    return n > 0.0f ? F3{a.x / n, a.y / n, a.z / n} : a;
}

F3 eyePosition(const Camera& c) {
    const float cp = std::cos(c.pitch);
    return {c.target[0] + c.distance * cp * std::sin(c.yaw),
            c.target[1] + c.distance * std::sin(c.pitch),
            c.target[2] + c.distance * cp * std::cos(c.yaw)};
}

// Base orthonormée de la caméra : f = direction du regard, s = droite, u = haut.
void cameraBasis(const Camera& c, F3& f, F3& s, F3& u) {
    const F3 eye = eyePosition(c);
    f = normalize(F3{c.target[0], c.target[1], c.target[2]} - eye);
    s = normalize(cross(f, F3{0.0f, 1.0f, 0.0f}));
    u = cross(s, f);
}

}  // namespace

void Camera::orbit(float dxPixels, float dyPixels) {
    yaw -= dxPixels * 0.005f;
    pitch = std::clamp(pitch + dyPixels * 0.005f, -1.5f, 1.5f);
}

void Camera::pan(float dxPixels, float dyPixels) {
    F3 f, s, u;
    cameraBasis(*this, f, s, u);
    const float k = distance * 0.0015f;
    target[0] += (-s.x * dxPixels + u.x * dyPixels) * k;
    target[1] += (-s.y * dxPixels + u.y * dyPixels) * k;
    target[2] += (-s.z * dxPixels + u.z * dyPixels) * k;
}

void Camera::zoom(float wheelTicks) {
    distance = std::clamp(distance * std::pow(0.9f, wheelTicks), 1.0f, 500.0f);
}

void Camera::viewProjection(float aspect, float out[16]) const {
    F3 f, s, u;
    cameraBasis(*this, f, s, u);
    const F3 eye = eyePosition(*this);

    // Vue (lookAt), stockage par colonnes.
    const float v[16] = {s.x,  u.x,  -f.x, 0.0f,
                         s.y,  u.y,  -f.y, 0.0f,
                         s.z,  u.z,  -f.z, 0.0f,
                         -dot(s, eye), -dot(u, eye), dot(f, eye), 1.0f};

    // Projection perspective.
    const float zNear = 0.1f, zFar = 1000.0f;
    const float t = 1.0f / std::tan(0.5f * fovY);
    const float p[16] = {t / aspect, 0.0f, 0.0f, 0.0f,
                         0.0f,       t,    0.0f, 0.0f,
                         0.0f,       0.0f, (zFar + zNear) / (zNear - zFar), -1.0f,
                         0.0f,       0.0f, 2.0f * zFar * zNear / (zNear - zFar), 0.0f};

    // out = p * v  (colonne c, ligne r -> indice c*4 + r)
    for (int c = 0; c < 4; ++c)
        for (int r = 0; r < 4; ++r) {
            float sum = 0.0f;
            for (int k = 0; k < 4; ++k) sum += p[k * 4 + r] * v[c * 4 + k];
            out[c * 4 + r] = sum;
        }
}

}  // namespace pl
