#include "Particle.h"

#include <cmath>

double coulomb_core(const Particle2D &a, const Particle2D &b) {
    const auto dx = a.x - b.x;
    const auto dy = a.y - b.y;

    const auto distance = std::sqrt(dx * dx + dy * dy);
    return a.q * b.q / distance + std::pow(distance, -8.0);
}

double coulomb_core(const Particle3D &a, const Particle3D &b) {
    const auto dx = a.x - b.x;
    const auto dy = a.y - b.y;
    const auto dz = a.z - b.z;

    const auto distance = std::sqrt(dx * dx + dy * dy + dz * dz);
    return a.q * b.q / distance + std::pow(distance, -8.0);
}
