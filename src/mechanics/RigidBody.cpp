#include "physicslab/mechanics/RigidBody.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>

#include "physicslab/core/Constants.hpp"
#include "physicslab/mechanics/Pendulum.hpp"   // ellipticK

namespace pl {

namespace inertia {

Vec3 sphere(double m, double r) {
    const double i = 0.4 * m * r * r;
    return {i, i, i};
}

Vec3 box(double m, double a, double b, double c) {
    const double k = m / 12.0;
    return {k * (b * b + c * c), k * (a * a + c * c), k * (a * a + b * b)};
}

Vec3 cylinder(double m, double r, double h) {
    const double transverse = m * (3.0 * r * r + h * h) / 12.0;
    return {transverse, transverse, 0.5 * m * r * r};
}

}  // namespace inertia

void jacobiSnCnDn(double u, double m, double& sn, double& cn, double& dn) {
    assert(m >= 0.0 && m < 1.0 && "jacobiSnCnDn : il faut 0 <= m < 1");
    if (m < 1e-14) {   // limite trigonométrique (le préfacteur 1/k de la série n'est plus défini)
        sn = std::sin(u);
        cn = std::cos(u);
        dn = 1.0;
        return;
    }
    // K = pi / (2 AGM(1, k')) et K' = pi / (2 AGM(1, k)), avec k = sqrt(m) et le module complémentaire k' = sqrt(1 - m) formé SANS passer par
    // 1 - k^2 refait deux fois (ellipticK(sqrt(1 - m)) perd la précision pour m petit : mesuré, nome faux de 0,5 % à m = 1e-14).
    const double k = std::sqrt(m), kc = std::sqrt(1.0 - m);
    auto agm = [](double a, double b) {
        for (int i = 0; i < 40 && std::abs(a - b) > 1e-16 * a; ++i) {
            const double next = 0.5 * (a + b);
            b = std::sqrt(a * b);
            a = next;
        }
        return a;
    };
    const double kk = 0.5 * constants::pi / agm(1.0, kc), kkPrime = 0.5 * constants::pi / agm(1.0, k);
    const double q = std::exp(-constants::pi * kkPrime / kk);
    const double zeta = constants::pi * u / (2.0 * kk);
    // Série géométrique de raison q : 41 / |ln q| termes ramènent q^n sous 1e-18.
    const int terms = std::clamp(static_cast<int>(std::ceil(41.0 / -std::log(q))) + 2, 4, 5000);
    double sSn = 0.0, sCn = 0.0, sDn = 0.0, qn = std::sqrt(q);   // qn = q^(n + 1/2)
    for (int n = 0; n < terms; ++n) {
        const double odd = (2 * n + 1) * zeta, q2 = qn * qn;       // q2 = q^(2n+1)
        sSn += qn / (1.0 - q2) * std::sin(odd);
        sCn += qn / (1.0 + q2) * std::cos(odd);
        qn *= q;
    }
    double qe = q;   // q^n, n >= 1
    for (int n = 1; n < terms; ++n) {
        sDn += qe / (1.0 + qe * qe) * std::cos(2.0 * n * zeta);
        qe *= q;
    }
    const double c = 2.0 * constants::pi / kk;
    sn = c / k * sSn;
    cn = c / k * sCn;
    dn = constants::pi / (2.0 * kk) + c * sDn;
}

State packRotation(const RotationState& s) { return {s.q.w, s.q.x, s.q.y, s.q.z, s.omega.x, s.omega.y, s.omega.z}; }

RotationState unpackRotation(const State& y) { return {Quaternion{y[0], y[1], y[2], y[3]}, Vec3{y[4], y[5], y[6]}}; }

Vec3 rotationOmegaDot(const Vec3& inertia, const Vec3& omega, const Vec3& torque) {
    return {((inertia.y - inertia.z) * omega.y * omega.z + torque.x) / inertia.x,
            ((inertia.z - inertia.x) * omega.z * omega.x + torque.y) / inertia.y,
            ((inertia.x - inertia.y) * omega.x * omega.y + torque.z) / inertia.z};
}

OdeFunction rotationRhs(const Vec3& inertia, const TorqueFunction& torque) {
    return [inertia, torque](double, const State& y, State& dydt) {
        const RotationState s = unpackRotation(y);
        // q' = 1/2 q (0, w)
        const Quaternion dq = s.q * Quaternion{0.0, s.omega.x, s.omega.y, s.omega.z};
        dydt[0] = 0.5 * dq.w;
        dydt[1] = 0.5 * dq.x;
        dydt[2] = 0.5 * dq.y;
        dydt[3] = 0.5 * dq.z;
        const Vec3 dw = rotationOmegaDot(inertia, s.omega, torque ? torque(s.q) : Vec3{});
        dydt[4] = dw.x;
        dydt[5] = dw.y;
        dydt[6] = dw.z;
    };
}

double FreeBodyProblem::energy(const RotationState& s) const {
    return 0.5 * (inertia.x * s.omega.x * s.omega.x + inertia.y * s.omega.y * s.omega.y + inertia.z * s.omega.z * s.omega.z);
}

Vec3 FreeBodyProblem::angularMomentumBody(const RotationState& s) const {
    return {inertia.x * s.omega.x, inertia.y * s.omega.y, inertia.z * s.omega.z};
}

Vec3 FreeBodyProblem::angularMomentumSpace(const RotationState& s) const { return s.q.rotate(angularMomentumBody(s)); }

RotationState FreeBodyProblem::exactSymmetric(double t) const {
    assert(std::abs(inertia.x - inertia.y) <= 1e-12 * inertia.x && "exactSymmetric : il faut I1 = I2");
    const Vec3 lBody = angularMomentumBody({q0, omega0});
    const double l = lBody.norm();
    const double omegaRel = (inertia.z - inertia.x) * omega0.z / inertia.x;   // vitesse de (w1, w2) dans le repère du corps
    const double c = std::cos(omegaRel * t), s = std::sin(omegaRel * t);
    RotationState out;
    out.omega = {omega0.x * c - omega0.y * s, omega0.x * s + omega0.y * c, omega0.z};
    const Quaternion spin = Quaternion::fromAxisAngle({0.0, 0.0, 1.0}, -omegaRel * t);
    if (l == 0.0) {
        out.q = q0 * spin;   // pas de rotation du tout (w = 0) : l'orientation ne change pas
    } else {
        const Quaternion precession = Quaternion::fromAxisAngle(q0.rotate(lBody), l / inertia.x * t);
        out.q = precession * q0 * spin;
    }
    return out;
}

namespace {

struct AsymmetricParams {
    bool nearAxis3 = true;
    double lambda = 0.0, m = 0.0, b = 0.0;
};

AsymmetricParams asymmetricParams(const Vec3& inertia, const Vec3& w0) {
    assert(inertia.x < inertia.y && inertia.y < inertia.z && "solution asymétrique : il faut I1 < I2 < I3");
    assert(w0.x > 0.0 && w0.y == 0.0 && w0.z > 0.0 && "solution asymétrique : départ (a, 0, c) avec a, c > 0");
    const double i1 = inertia.x, i2 = inertia.y, i3 = inertia.z;
    const double twoE = i1 * w0.x * w0.x + i3 * w0.z * w0.z;
    const double l2 = i1 * i1 * w0.x * w0.x + i3 * i3 * w0.z * w0.z;
    AsymmetricParams p;
    const double mLandau = (i2 - i1) * (twoE * i3 - l2) / ((i3 - i2) * (l2 - twoE * i1));
    p.nearAxis3 = twoE * i2 < l2;
    if (p.nearAxis3) {
        p.lambda = w0.z * std::sqrt((i3 - i2) * (i3 - i1) / (i1 * i2));
        p.m = mLandau;
    } else {
        p.lambda = w0.x * std::sqrt((i3 - i1) * (i2 - i1) / (i2 * i3));
        p.m = 1.0 / mLandau;
    }
    p.b = (i3 - i1) * w0.x * w0.z / (i2 * p.lambda);
    return p;
}

}  // namespace

Vec3 FreeBodyProblem::exactAsymmetricOmega(double t) const {
    const AsymmetricParams p = asymmetricParams(inertia, omega0);
    double sn, cn, dn;
    jacobiSnCnDn(p.lambda * t, p.m, sn, cn, dn);
    return p.nearAxis3 ? Vec3{omega0.x * cn, p.b * sn, omega0.z * dn} : Vec3{omega0.x * dn, p.b * sn, omega0.z * cn};
}

double FreeBodyProblem::asymmetricParameter() const { return asymmetricParams(inertia, omega0).m; }

double FreeBodyProblem::asymmetricPeriod() const {
    const AsymmetricParams p = asymmetricParams(inertia, omega0);
    return 4.0 * ellipticK(std::sqrt(p.m)) / p.lambda;
}

double FreeBodyProblem::intermediateAxisGrowthRate() const {
    return std::abs(omega0.y) * std::sqrt((inertia.z - inertia.y) * (inertia.y - inertia.x) / (inertia.x * inertia.z));
}

RotationState FreeBodyProblem::reference(double t) const {
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    State y = initialState();
    advance(rk, rhs(), 0.0, y, t, 10000000);
    return unpackRotation(y);
}

Quaternion quaternionExp(const Vec3& v) { return Quaternion::fromAxisAngle(v, v.norm()); }

void EulerRotation::step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) {
    const Vec3 a = rotationOmegaDot(inertia, s.omega, torque ? torque(s.q) : Vec3{});
    const Quaternion dq = s.q * Quaternion{0.0, s.omega.x, s.omega.y, s.omega.z};
    s.q = {s.q.w + 0.5 * h * dq.w, s.q.x + 0.5 * h * dq.x, s.q.y + 0.5 * h * dq.y, s.q.z + 0.5 * h * dq.z};
    s.omega += h * a;
}

void RK4Rotation::step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) {
    State y = packRotation(s);
    rk4_.step(rotationRhs(inertia, torque), 0.0, y, h);
    s = unpackRotation(y);
    if (renormalize_) s.q = s.q.normalized();
}

void LieHeunRotation::step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) {
    const Vec3 a1 = rotationOmegaDot(inertia, s.omega, torque ? torque(s.q) : Vec3{});
    // Prédicteur : orientation et vitesse d'Euler (par l'exponentielle).
    const Quaternion qPredicted = s.q * quaternionExp(h * s.omega);
    const Vec3 wPredicted = s.omega + h * a1;
    const Vec3 a2 = rotationOmegaDot(inertia, wPredicted, torque ? torque(qPredicted) : Vec3{});
    // Correcteur : la vitesse moyenne du pas (règle des trapèzes) pilote la rotation.
    const Vec3 wNew = s.omega + 0.5 * h * (a1 + a2);
    s.q = (s.q * quaternionExp(0.5 * h * (s.omega + wNew))).normalized();
    s.omega = wNew;
}

void SplittingRotation::step(const Vec3& inertia, const TorqueFunction& torque, RotationState& s, double h) {
    double pi[3] = {inertia.x * s.omega.x, inertia.y * s.omega.y, inertia.z * s.omega.z};
    const double moment[3] = {inertia.x, inertia.y, inertia.z};
    Quaternion q = s.q;

    auto kick = [&](double dt) {
        if (!torque) return;
        const Vec3 tau = torque(q);
        pi[0] += dt * tau.x;
        pi[1] += dt * tau.y;
        pi[2] += dt * tau.z;
    };
    // Rotateur d'axe k pendant tau : w_k constant ; (pi_a, pi_b) tourne de -theta, q est multiplié par la rotation de +theta autour de e_k.
    auto rotor = [&](int k, double tau) {
        const int a = (k + 1) % 3, b = (k + 2) % 3;
        const double theta = pi[k] / moment[k] * tau, c = std::cos(theta), sn = std::sin(theta);
        const double pa = pi[a], pb = pi[b];
        pi[a] = pa * c + pb * sn;
        pi[b] = -pa * sn + pb * c;
        const Vec3 axis = k == 0 ? Vec3{1.0, 0.0, 0.0} : (k == 1 ? Vec3{0.0, 1.0, 0.0} : Vec3{0.0, 0.0, 1.0});
        q = q * Quaternion::fromAxisAngle(axis, theta);
    };

    kick(0.5 * h);
    rotor(0, 0.5 * h);
    rotor(1, 0.5 * h);
    rotor(2, h);
    rotor(1, 0.5 * h);
    rotor(0, 0.5 * h);
    kick(0.5 * h);

    s.q = q.normalized();   // ne retire que l'accumulation des arrondis (de l'ordre de 1e-16 par pas)
    s.omega = {pi[0] / moment[0], pi[1] / moment[1], pi[2] / moment[2]};
}

RotationState integrateRotation(RotationIntegrator& integrator, const Vec3& inertia, const TorqueFunction& torque, RotationState s, double tEnd,
                                int steps) {
    const double h = tEnd / steps;
    for (int i = 0; i < steps; ++i) integrator.step(inertia, torque, s, h);
    return s;
}

double rotationDistance(const RotationState& a, const RotationState& b) {
    const Mat3 ra = a.q.normalized().toMat3(), rb = b.q.normalized().toMat3();
    double sum = (a.omega - b.omega).norm2();
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) sum += (ra.m[i][j] - rb.m[i][j]) * (ra.m[i][j] - rb.m[i][j]);
    return std::sqrt(sum);
}

TorqueFunction HeavyTopProblem::torque() const {
    const double m = mass, l = lever, g = gravity;
    return [m, l, g](const Quaternion& q) {
        const Vec3 gBody = q.conjugate().rotate({0.0, 0.0, -g});   // pesanteur vue dans le repère du corps
        return m * l * cross(Vec3{0.0, 0.0, 1.0}, gBody);
    };
}

double HeavyTopProblem::energy(const RotationState& s) const {
    return 0.5 * (inertia.x * s.omega.x * s.omega.x + inertia.y * s.omega.y * s.omega.y + inertia.z * s.omega.z * s.omega.z) +
           mass * gravity * lever * cosTheta(s.q);
}

Vec3 HeavyTopProblem::angularMomentumSpace(const RotationState& s) const {
    return s.q.rotate({inertia.x * s.omega.x, inertia.y * s.omega.y, inertia.z * s.omega.z});
}

RotationState HeavyTopProblem::reference(double t) const {
    RK45 rk;
    rk.relTol = 1e-13;
    rk.absTol = 1e-15;
    State y = initialState();
    advance(rk, rhs(), 0.0, y, t, 10000000);
    return unpackRotation(y);
}

SteadyPrecession HeavyTopProblem::steadyPrecession(double theta, double spin, bool slow) const {
    SteadyPrecession sp;
    sp.theta = theta;
    sp.spin = spin;
    const double u = std::cos(theta), a = inertia.z * spin, mgl = mass * gravity * lever;
    const double discriminant = a * a - 4.0 * inertia.x * mgl * u;
    if (spin <= 0.0 || discriminant < 0.0) return sp;
    const double root = std::sqrt(discriminant);
    if (slow) {
        sp.phiDot = 2.0 * mgl / (a + root);          // (a - root) / (2 I1 u) sous une forme stable, valable aussi pour u = 0
    } else {
        if (std::abs(u) < 1e-12) return sp;          // cos theta = 0 : la branche rapide part à l'infini
        sp.phiDot = (a + root) / (2.0 * inertia.x * u);
    }
    sp.psiDot = spin - sp.phiDot * u;
    sp.valid = true;
    return sp;
}

void HeavyTopProblem::startSteady(const SteadyPrecession& sp) {
    q0 = Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, sp.theta);
    omega0 = {0.0, sp.phiDot * std::sin(sp.theta), sp.spin};
}

RotationState HeavyTopProblem::exactSteady(const SteadyPrecession& sp, double t) const {
    RotationState s;
    s.q = Quaternion::fromAxisAngle({0.0, 0.0, 1.0}, sp.phiDot * t) * Quaternion::fromAxisAngle({1.0, 0.0, 0.0}, sp.theta) *
          Quaternion::fromAxisAngle({0.0, 0.0, 1.0}, sp.psiDot * t);
    const double psi = sp.psiDot * t, sinT = std::sin(sp.theta);
    s.omega = {sp.phiDot * sinT * std::sin(psi), sp.phiDot * sinT * std::cos(psi), sp.spin};
    return s;
}

NutationRange HeavyTopProblem::nutation() const {
    NutationRange out;
    const RotationState s0{q0, omega0};
    const double i1 = inertia.x, mgl = mass * gravity * lever;
    const double l3 = axialMomentum(s0), lz = verticalMomentum(s0);
    const double ePrime = energy(s0) - l3 * l3 / (2.0 * inertia.z);
    // f(u) = c3 u^3 + c2 u^2 + c1 u + c0
    const double c3 = 2.0 * mgl / i1, c2 = -(2.0 * ePrime / i1 + l3 * l3 / (i1 * i1));
    const double c1 = -2.0 * mgl / i1 + 2.0 * lz * l3 / (i1 * i1), c0 = 2.0 * ePrime / i1 - lz * lz / (i1 * i1);
    if (c3 <= 0.0) return out;
    // Racines du cubique unitaire u^3 + A u^2 + B u + C par la méthode trigonométrique (trois racines réelles), affinées par Newton.
    const double A = c2 / c3, B = c1 / c3, C = c0 / c3;
    const double q = (A * A - 3.0 * B) / 9.0, r = (2.0 * A * A * A - 9.0 * A * B + 27.0 * C) / 54.0;
    if (q <= 0.0 || r * r >= q * q * q) return out;
    const double phi = std::acos(r / std::sqrt(q * q * q));
    double roots[3];
    for (int k = 0; k < 3; ++k) {
        double u = -2.0 * std::sqrt(q) * std::cos((phi + 2.0 * constants::pi * k) / 3.0) - A / 3.0;
        for (int it = 0; it < 4; ++it) {
            const double f = ((u + A) * u + B) * u + C, df = (3.0 * u + 2.0 * A) * u + B;
            if (df != 0.0) u -= f / df;
        }
        roots[k] = u;
    }
    std::sort(roots, roots + 3);
    out.u1 = roots[0];
    out.u2 = roots[1];
    out.u3 = roots[2];
    out.valid = out.u1 < out.u2 && out.u2 < out.u3;
    return out;
}

double HeavyTopProblem::nutationPeriod() const {
    const NutationRange n = nutation();
    if (!n.valid) return 0.0;
    const double beta = 2.0 * mass * gravity * lever / inertia.x;
    return 4.0 * ellipticK(std::sqrt((n.u2 - n.u1) / (n.u3 - n.u1))) / std::sqrt(beta * (n.u3 - n.u1));
}

double HeavyTopProblem::sleepingCriticalSpin() const { return 2.0 * std::sqrt(inertia.x * mass * gravity * lever) / inertia.z; }

}  // namespace pl
