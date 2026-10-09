#include "physicslab/render/FieldPlot.hpp"

#include <algorithm>
#include <cmath>

namespace pl {

void divergingColor(double t, float rgb[3]) {
    const float a = static_cast<float>(std::clamp(std::abs(t), 0.0, 1.0));
    static const float zero[3] = {0.42f, 0.45f, 0.50f};
    static const float negative[3] = {0.20f, 0.45f, 1.00f};
    static const float positive[3] = {1.00f, 0.32f, 0.22f};
    const float* end = t < 0.0 ? negative : positive;
    for (int c = 0; c < 3; ++c) rgb[c] = zero[c] + (end[c] - zero[c]) * a;
}

std::vector<Vertex> makeProfile(const std::vector<double>& u, double dx, float heightScale, float colorScale, float z) {
    std::vector<Vertex> v;
    const int n = static_cast<int>(u.size());
    v.reserve(n);
    const double half = 0.5 * (n - 1);
    for (int i = 0; i < n; ++i) {
        float rgb[3];
        divergingColor(u[i] / colorScale, rgb);
        v.push_back({static_cast<float>((i - half) * dx), heightScale * static_cast<float>(u[i]), z, rgb[0], rgb[1], rgb[2]});
    }
    return v;
}

std::vector<Vertex> makeRelief(const std::vector<double>& u, int nx, int ny, double dx, float heightScale, float colorScale, int stride) {
    stride = std::max(stride, 1);
    const double hx = 0.5 * (nx - 1), hz = 0.5 * (ny - 1);
    auto vertex = [&](int i, int j) {
        const double value = u[static_cast<std::size_t>(j) * nx + i];
        float rgb[3];
        divergingColor(value / colorScale, rgb);
        return Vertex{static_cast<float>((i - hx) * dx), heightScale * static_cast<float>(value),
                      static_cast<float>((j - hz) * dx), rgb[0], rgb[1], rgb[2]};
    };

    std::vector<Vertex> v;
    // lignes parallèles à x (une rangée sur `stride`), puis parallèles à z ; les deux dernières lignes du bord sont toujours tracées
    for (int j = 0; j < ny; ++j) {
        if (j % stride != 0 && j != ny - 1) continue;
        for (int i = 0; i + 1 < nx; ++i) {
            v.push_back(vertex(i, j));
            v.push_back(vertex(i + 1, j));
        }
    }
    for (int i = 0; i < nx; ++i) {
        if (i % stride != 0 && i != nx - 1) continue;
        for (int j = 0; j + 1 < ny; ++j) {
            v.push_back(vertex(i, j));
            v.push_back(vertex(i, j + 1));
        }
    }
    return v;
}

}  // namespace pl
