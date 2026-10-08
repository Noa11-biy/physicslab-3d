// Vecteur 3D en double précision : brique de base de toute la physique du moteur.
#pragma once

#include <cmath>

namespace pl {

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    constexpr Vec3() = default;
    constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    constexpr Vec3 operator-() const { return {-x, -y, -z}; }
    constexpr Vec3& operator+=(const Vec3& o) { x += o.x; y += o.y; z += o.z; return *this; }
    constexpr Vec3& operator-=(const Vec3& o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    constexpr Vec3& operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
    constexpr Vec3& operator/=(double s) { x /= s; y /= s; z /= s; return *this; }

    constexpr double norm2() const { return x * x + y * y + z * z; }
    double norm() const { return std::sqrt(norm2()); }

    // Vecteur unitaire ; le vecteur nul reste nul (pas de division par zéro).
    Vec3 normalized() const {
        const double n = norm();
        return n > 0.0 ? Vec3{x / n, y / n, z / n} : Vec3{};
    }
};

constexpr Vec3 operator+(Vec3 a, const Vec3& b) { return a += b; }
constexpr Vec3 operator-(Vec3 a, const Vec3& b) { return a -= b; }
constexpr Vec3 operator*(Vec3 a, double s) { return a *= s; }
constexpr Vec3 operator*(double s, Vec3 a) { return a *= s; }
constexpr Vec3 operator/(Vec3 a, double s) { return a /= s; }

constexpr double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

constexpr Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

}  // namespace pl
