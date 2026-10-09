// Matrice 3x3 en double (tenseur d'inertie, rotations, tenseurs d'ordre 2 de base).
// Stockage ligne par ligne : m[ligne][colonne].
#pragma once

#include <cassert>
#include <cmath>

#include "physicslab/core/Vec3.hpp"

namespace pl {

struct Mat3 {
    double m[3][3] = {};  // matrice nulle par défaut

    static constexpr Mat3 identity() { return diagonal(1.0, 1.0, 1.0); }

    static constexpr Mat3 diagonal(double a, double b, double c) {
        Mat3 r;
        r.m[0][0] = a;
        r.m[1][1] = b;
        r.m[2][2] = c;
        return r;
    }

    // Matrice antisymétrique [v]x telle que [v]x * w = v x w.
    static constexpr Mat3 skew(const Vec3& v) {
        Mat3 r;
        r.m[0][1] = -v.z; r.m[0][2] = v.y;
        r.m[1][0] = v.z;  r.m[1][2] = -v.x;
        r.m[2][0] = -v.y; r.m[2][1] = v.x;
        return r;
    }

    constexpr Mat3 transposed() const {
        Mat3 r;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) r.m[i][j] = m[j][i];
        return r;
    }

    constexpr double determinant() const {
        return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    }

    // Inverse par la comatrice. Précondition : déterminant non nul.
    Mat3 inverse() const {
        const double d = determinant();
        assert(std::abs(d) > 0.0 && "Mat3::inverse sur une matrice singulière");
        const double id = 1.0 / d;
        Mat3 r;
        r.m[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * id;
        r.m[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id;
        r.m[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
        r.m[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id;
        r.m[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
        r.m[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * id;
        r.m[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * id;
        r.m[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * id;
        r.m[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * id;
        return r;
    }
};

constexpr Vec3 operator*(const Mat3& a, const Vec3& v) {
    return {a.m[0][0] * v.x + a.m[0][1] * v.y + a.m[0][2] * v.z,
            a.m[1][0] * v.x + a.m[1][1] * v.y + a.m[1][2] * v.z,
            a.m[2][0] * v.x + a.m[2][1] * v.y + a.m[2][2] * v.z};
}

constexpr Mat3 operator*(const Mat3& a, const Mat3& b) {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) r.m[i][j] += a.m[i][k] * b.m[k][j];
    return r;
}

constexpr Mat3 operator*(double s, const Mat3& a) {
    Mat3 r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) r.m[i][j] = s * a.m[i][j];
    return r;
}

}  // namespace pl
