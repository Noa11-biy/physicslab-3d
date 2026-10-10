#include "physicslab/waves/Wave.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

namespace pl::waves {

double spongeSigmaMax(double speed, double depth, double reflection, int order) {
    assert(speed > 0.0 && depth > 0.0 && reflection > 0.0 && reflection < 1.0 && order >= 1);
    return (order + 1) * speed * std::log(1.0 / reflection) / (2.0 * depth);
}

namespace {

// Profil de l'éponge : sigma_max s^order, s = 1 contre le mur, 0 à l'entrée de la couche. Les premières dérivées de sigma sont
// nulles à l'entrée jusqu'à l'ordre order - 1 : plus l'ordre est grand, moins l'entrée réfléchit.
double spongeProfile(int distanceFromWall, int layer, double sigmaMax, int order) {
    if (distanceFromWall >= layer) return 0.0;
    const double s = static_cast<double>(layer - distanceFromWall) / static_cast<double>(layer);
    return sigmaMax * std::pow(s, order);
}

}  // namespace

// ------------------------------------------------------------------ 1D ------

Wave1D::Wave1D(const Wave1DParams& params) : p_(params) {
    assert(p_.cells >= 2 && p_.length > 0.0 && p_.speed > 0.0 && p_.cfl > 0.0);
    p_.spongeCells = std::clamp(p_.spongeCells, 1, p_.cells / 2);
    dx_ = p_.length / p_.cells;
    dt_ = p_.cfl * dx_ / p_.speed;

    const std::size_t n = static_cast<std::size_t>(p_.cells) + 1;
    prev_.assign(n, 0.0);
    cur_.assign(n, 0.0);
    next_.assign(n, 0.0);
    sigma_.assign(n, 0.0);

    const double sigmaMax = spongeSigmaMax(p_.speed, p_.spongeCells * dx_, p_.spongeReflection, p_.spongeOrder);
    for (int i = 0; i <= p_.cells; ++i) {
        if (p_.left == Edge::Absorbing) sigma_[i] += spongeProfile(i, p_.spongeCells, sigmaMax, p_.spongeOrder);
        if (p_.right == Edge::Absorbing) sigma_[i] += spongeProfile(p_.cells - i, p_.spongeCells, sigmaMax, p_.spongeOrder);
    }
}

void Wave1D::setInitial(const std::function<double(double)>& u0, const std::function<double(double)>& v0) {
    const int n = p_.cells;
    const double c2 = p_.cfl * p_.cfl;
    const bool leftFixed = p_.left != Edge::Free, rightFixed = p_.right != Edge::Free;

    std::vector<double> v(n + 1, 0.0);
    for (int i = 0; i <= n; ++i) {
        cur_[i] = u0 ? u0(x(i)) : 0.0;
        if (v0) v[i] = v0(x(i));
    }
    if (leftFixed) { cur_[0] = 0.0; v[0] = 0.0; }
    if (rightFixed) { cur_[n] = 0.0; v[n] = 0.0; }

    // u^{-1} = u^0 - dt v0 + ½ dt² a0, avec a0 = c² u_xx - 2 sigma v0 (Taylor d'ordre 2 en arrière)
    for (int i = 0; i <= n; ++i) {
        double lap;
        if (i == 0) lap = leftFixed ? 0.0 : 2.0 * (cur_[1] - cur_[0]);
        else if (i == n) lap = rightFixed ? 0.0 : 2.0 * (cur_[n - 1] - cur_[n]);
        else lap = cur_[i - 1] - 2.0 * cur_[i] + cur_[i + 1];
        prev_[i] = cur_[i] - dt_ * v[i] + 0.5 * c2 * lap - sigma_[i] * dt_ * dt_ * v[i];
    }
    if (leftFixed) prev_[0] = 0.0;
    if (rightFixed) prev_[n] = 0.0;
    steps_ = 0;
}

void Wave1D::setStates(const std::function<double(double)>& previous, const std::function<double(double)>& current) {
    const int n = p_.cells;
    for (int i = 0; i <= n; ++i) {
        prev_[i] = previous(x(i));
        cur_[i] = current(x(i));
    }
    if (p_.left != Edge::Free) { prev_[0] = 0.0; cur_[0] = 0.0; }
    if (p_.right != Edge::Free) { prev_[n] = 0.0; cur_[n] = 0.0; }
    steps_ = 0;
}

void Wave1D::step() {
    const int n = p_.cells;
    const double c2 = p_.cfl * p_.cfl;
    for (int i = 0; i <= n; ++i) {
        double lap;
        if (i == 0) {
            if (p_.left != Edge::Free) { next_[0] = 0.0; continue; }
            lap = 2.0 * (cur_[1] - cur_[0]);  // point fantôme u_{-1} = u_1 : dérivée nulle
        } else if (i == n) {
            if (p_.right != Edge::Free) { next_[n] = 0.0; continue; }
            lap = 2.0 * (cur_[n - 1] - cur_[n]);
        } else {
            lap = cur_[i - 1] - 2.0 * cur_[i] + cur_[i + 1];
        }
        const double s = sigma_[i] * dt_;
        next_[i] = (2.0 * cur_[i] - (1.0 - s) * prev_[i] + c2 * lap) / (1.0 + s);
    }
    prev_.swap(next_);  // prev = u^{n+1}, next = ancien u^{n-1} (tampon libre)
    prev_.swap(cur_);   // prev = u^n, cur = u^{n+1}
    ++steps_;
}

double Wave1D::energy() const {
    const int n = p_.cells;
    double kinetic = 0.0;
    for (int i = 0; i <= n; ++i) {
        const double w = (i == 0 || i == n) ? 0.5 : 1.0;
        const double v = (cur_[i] - prev_[i]) / dt_;
        kinetic += w * dx_ * v * v;
    }
    double potential = 0.0;
    for (int j = 0; j < n; ++j) potential += (cur_[j + 1] - cur_[j]) * (prev_[j + 1] - prev_[j]) / dx_;
    return 0.5 * kinetic + 0.5 * p_.speed * p_.speed * potential;
}

double Wave1D::maxAbs() const {
    double m = 0.0;
    for (double v : cur_) {
        const double a = std::abs(v);
        if (!(a <= m)) m = a;  // contrairement à std::max, un NaN se propage : un calcul qui a explosé ne passe pas pour un calcul calme
    }
    return m;
}

// ------------------------------------------------------------------ 2D ------

Wave2D::Wave2D(const Wave2DParams& params) : p_(params) {
    assert(p_.cellsX >= 2 && p_.cellsY >= 2 && p_.dx > 0.0 && p_.speed > 0.0 && p_.cfl > 0.0);
    p_.spongeCells = std::clamp(p_.spongeCells, 1, std::min(p_.cellsX, p_.cellsY) / 2);

    const std::size_t count = static_cast<std::size_t>(pointsX()) * static_cast<std::size_t>(pointsY());
    prev_.assign(count, 0.0);
    cur_.assign(count, 0.0);
    next_.assign(count, 0.0);
    sigma_.assign(count, 0.0);

    // Vitesse en chaque noeud ; dt est fixé par la vitesse maximale.
    assert(p_.speedField.empty() || p_.speedField.size() == count);
    std::vector<double> speed(count, p_.speed);
    if (!p_.speedField.empty()) speed = p_.speedField;
    cMax_ = *std::max_element(speed.begin(), speed.end());
    assert(cMax_ > 0.0);
    dt_ = p_.cfl * p_.dx / cMax_;
    ratio_ = dt_ * dt_ / (p_.dx * p_.dx);

    // Coefficients des faces : (dt² / dx²) fois la moyenne arithmétique des c² des deux noeuds (c = vitesse uniforme : C²).
    kx_.assign(count, 0.0);
    ky_.assign(count, 0.0);
    const int stride = p_.cellsX + 1;
    for (int j = 0; j <= p_.cellsY; ++j)
        for (int i = 0; i <= p_.cellsX; ++i) {
            const int k = index(i, j);
            const double c2 = speed[k] * speed[k];
            if (i < p_.cellsX) kx_[k] = ratio_ * 0.5 * (c2 + speed[k + 1] * speed[k + 1]);
            if (j < p_.cellsY) ky_[k] = ratio_ * 0.5 * (c2 + speed[k + stride] * speed[k + stride]);
        }

    // Amortissement additif : sigma(x, y) = sigma_x(x) + sigma_y(y). Aux coins les deux couches se cumulent.
    const double sigmaMax = spongeSigmaMax(cMax_, p_.spongeCells * p_.dx, p_.spongeReflection, p_.spongeOrder);
    const int order = p_.spongeOrder;
    for (int j = 0; j <= p_.cellsY; ++j) {
        for (int i = 0; i <= p_.cellsX; ++i) {
            double s = 0.0;
            if (p_.left == Edge::Absorbing) s += spongeProfile(i, p_.spongeCells, sigmaMax, order);
            if (p_.right == Edge::Absorbing) s += spongeProfile(p_.cellsX - i, p_.spongeCells, sigmaMax, order);
            if (p_.bottom == Edge::Absorbing) s += spongeProfile(j, p_.spongeCells, sigmaMax, order);
            if (p_.top == Edge::Absorbing) s += spongeProfile(p_.cellsY - j, p_.spongeCells, sigmaMax, order);
            sigma_[index(i, j)] = s;
        }
    }
}

std::vector<double> Wave2D::makeSpeedField(const Wave2DParams& g, const std::function<double(double, double)>& c) {
    std::vector<double> field(static_cast<std::size_t>(g.cellsX + 1) * static_cast<std::size_t>(g.cellsY + 1));
    for (int j = 0; j <= g.cellsY; ++j)
        for (int i = 0; i <= g.cellsX; ++i) field[static_cast<std::size_t>(j) * (g.cellsX + 1) + i] = c(i * g.dx, j * g.dx);
    return field;
}

double Wave2D::weightedLaplacian(const std::vector<double>& f, int i, int j) const {
    const int nx = p_.cellsX, ny = p_.cellsY, stride = nx + 1;
    const int k = index(i, j);
    const double c = f[k];
    // Flux de chaque face. Au-delà d'un bord libre : point fantôme symétrique, donc face symétrique (même coefficient, valeur opposée) ;
    // sur un bord fixe le noeud n'est jamais mis à jour (flux inutilisé).
    const double west = i > 0 ? kx_[k - 1] * (f[k - 1] - c) : (p_.left == Edge::Free ? kx_[k] * (f[k + 1] - c) : 0.0);
    const double east = i < nx ? kx_[k] * (f[k + 1] - c) : (p_.right == Edge::Free ? kx_[k - 1] * (f[k - 1] - c) : 0.0);
    const double south = j > 0 ? ky_[k - stride] * (f[k - stride] - c) : (p_.bottom == Edge::Free ? ky_[k] * (f[k + stride] - c) : 0.0);
    const double north = j < ny ? ky_[k] * (f[k + stride] - c) : (p_.top == Edge::Free ? ky_[k - stride] * (f[k - stride] - c) : 0.0);
    return west + east + south + north;
}

void Wave2D::clampBoundaries(std::vector<double>& f) const {
    const int nx = p_.cellsX, ny = p_.cellsY;
    for (int j = 0; j <= ny; ++j) {
        if (p_.left != Edge::Free) f[index(0, j)] = 0.0;
        if (p_.right != Edge::Free) f[index(nx, j)] = 0.0;
    }
    for (int i = 0; i <= nx; ++i) {
        if (p_.bottom != Edge::Free) f[index(i, 0)] = 0.0;
        if (p_.top != Edge::Free) f[index(i, ny)] = 0.0;
    }
}

void Wave2D::setInitial(const std::function<double(double, double)>& u0, const std::function<double(double, double)>& v0) {
    const int nx = p_.cellsX, ny = p_.cellsY;
    std::vector<double> v(cur_.size(), 0.0);
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const double x = i * p_.dx, y = j * p_.dx;
            cur_[index(i, j)] = u0 ? u0(x, y) : 0.0;
            if (v0) v[index(i, j)] = v0(x, y);
        }
    }
    clampBoundaries(cur_);
    clampBoundaries(v);

    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const int k = index(i, j);
            prev_[k] = cur_[k] - dt_ * v[k] + 0.5 * weightedLaplacian(cur_, i, j) - sigma_[k] * dt_ * dt_ * v[k];
        }
    }
    clampBoundaries(prev_);
    steps_ = 0;
}

void Wave2D::setStates(const std::function<double(double, double)>& previous, const std::function<double(double, double)>& current) {
    for (int j = 0; j <= p_.cellsY; ++j) {
        for (int i = 0; i <= p_.cellsX; ++i) {
            const double x = i * p_.dx, y = j * p_.dx;
            prev_[index(i, j)] = previous(x, y);
            cur_[index(i, j)] = current(x, y);
        }
    }
    clampBoundaries(prev_);
    clampBoundaries(cur_);
    steps_ = 0;
}

void Wave2D::step() {
    const int nx = p_.cellsX, ny = p_.cellsY;
    for (int j = 0; j <= ny; ++j) {
        const bool rowFixed = (j == 0 && p_.bottom != Edge::Free) || (j == ny && p_.top != Edge::Free);
        for (int i = 0; i <= nx; ++i) {
            const int k = index(i, j);
            if (rowFixed || (i == 0 && p_.left != Edge::Free) || (i == nx && p_.right != Edge::Free)) {
                next_[k] = 0.0;
                continue;
            }
            const double s = sigma_[k] * dt_;
            next_[k] = (2.0 * cur_[k] - (1.0 - s) * prev_[k] + weightedLaplacian(cur_, i, j)) / (1.0 + s);
        }
    }
    // sources : force A sin(w t) ramp(t) au noeud le plus proche (jamais sur un bord fixe)
    const double t = steps_ * dt_;
    for (const PointSource& src : sources_) {
        const int i = std::clamp(static_cast<int>(std::lround(src.x / p_.dx)), 0, nx), j = std::clamp(static_cast<int>(std::lround(src.y / p_.dx)), 0, ny);
        const bool fixedNode = (i == 0 && p_.left != Edge::Free) || (i == nx && p_.right != Edge::Free) || (j == 0 && p_.bottom != Edge::Free) ||
                               (j == ny && p_.top != Edge::Free);
        if (fixedNode) continue;
        const double period = 6.283185307179586 / src.omega;
        const double r = src.rampPeriods > 0.0 ? std::min(1.0, t / (src.rampPeriods * period)) : 1.0;
        const double ramp = std::sin(0.5 * 3.14159265358979323846 * r) * std::sin(0.5 * 3.14159265358979323846 * r);
        next_[index(i, j)] += ratio_ * src.amplitude * ramp * std::sin(src.omega * t);
    }
    prev_.swap(next_);
    prev_.swap(cur_);
    ++steps_;
}

double Wave2D::energy() const {
    const int nx = p_.cellsX, ny = p_.cellsY;
    const double dx = p_.dx;
    auto weight = [](int i, int n) { return (i == 0 || i == n) ? 0.5 : 1.0; };

    double kinetic = 0.0;
    for (int j = 0; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            const int k = index(i, j);
            const double v = (cur_[k] - prev_[k]) / dt_;
            kinetic += weight(i, nx) * weight(j, ny) * dx * dx * v * v;
        }
    }
    // Arêtes horizontales (largeur transverse dx w_j) et verticales (dx w_i) : produit des différences nouvelle × ancienne.
    // c²_face = coefficient du schéma / (dt² / dx²) : la même moyenne que dans le schéma, sans quoi l'énergie ne serait plus conservée exactement.
    double potential = 0.0;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const int k = index(i, j);
            potential += weight(j, ny) * (kx_[k] / ratio_) * (cur_[k + 1] - cur_[k]) * (prev_[k + 1] - prev_[k]);
        }
    const int stride = nx + 1;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            const int k = index(i, j);
            potential += weight(i, nx) * (ky_[k] / ratio_) * (cur_[k + stride] - cur_[k]) * (prev_[k + stride] - prev_[k]);
        }
    return 0.5 * kinetic + 0.5 * potential;
}

double Wave2D::maxAbs() const {
    double m = 0.0;
    for (double v : cur_) {
        const double a = std::abs(v);
        if (!(a <= m)) m = a;  // un NaN se propage (voir Wave1D::maxAbs)
    }
    return m;
}

}  // namespace pl::waves
