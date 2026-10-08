// Quaternion unitaire pour représenter les orientations (corps rigides, M6).
// Convention de Hamilton : q = w + x i + y j + z k.
#pragma once

#include <cmath>

#include "physicslab/core/Mat3.hpp"
#include "physicslab/core/Vec3.hpp"

namespace pl {

struct Quaternion {
    double w = 1.0, x = 0.0, y = 0.0, z = 0.0;

    static Quaternion fromAxisAngle(const Vec3& axis, double angle) {
        const Vec3 u = axis.normalized();
        const double s = std::sin(0.5 * angle);
        return {std::cos(0.5 * angle), u.x * s, u.y * s, u.z * s};
    }

    constexpr Quaternion conjugate() const { return {w, -x, -y, -z}; }
    double norm() const { return std::sqrt(w * w + x * x + y * y + z * z); }

    // À rappeler régulièrement lors de l'intégration pour contrer la dérive numérique.
    Quaternion normalized() const {
        const double n = norm();
        return n > 0.0 ? Quaternion{w / n, x / n, y / n, z / n} : Quaternion{};
    }

    // Rotation du vecteur v : q v q*  (forme développée, sans produit quaternionique complet).
    Vec3 rotate(const Vec3& v) const {
        const Vec3 q{x, y, z};
        const Vec3 t = 2.0 * cross(q, v);
        return v + w * t + cross(q, t);
    }

    // Matrice de rotation équivalente (quaternion supposé unitaire).
    Mat3 toMat3() const {
        const double xx = x * x, yy = y * y, zz = z * z;
        const double xy = x * y, xz = x * z, yz = y * z;
        const double wx = w * x, wy = w * y, wz = w * z;
        Mat3 r;
        r.m[0][0] = 1.0 - 2.0 * (yy + zz); r.m[0][1] = 2.0 * (xy - wz);       r.m[0][2] = 2.0 * (xz + wy);
        r.m[1][0] = 2.0 * (xy + wz);       r.m[1][1] = 1.0 - 2.0 * (xx + zz); r.m[1][2] = 2.0 * (yz - wx);
        r.m[2][0] = 2.0 * (xz - wy);       r.m[2][1] = 2.0 * (yz + wx);       r.m[2][2] = 1.0 - 2.0 * (xx + yy);
        return r;
    }
};

// Produit de Hamilton : (a * b) applique d'abord b, puis a.
constexpr Quaternion operator*(const Quaternion& a, const Quaternion& b) {
    return {a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
            a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
            a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
            a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

}  // namespace pl
