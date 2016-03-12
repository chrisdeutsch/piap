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

double avg_pair_dist(const CanonicalEnsemble<Particle2D>::State &state) {
    const auto N = state.size();
    double distance = 0.0;
    for (auto it = state.cbegin(), end = state.cend(); it != end; ++it) {
        for (auto it2 = state.cbegin(); it2 != it; ++it2) {
            const auto dx = it->x - it2->x;
            const auto dy = it->y - it2->y;
            distance += std::sqrt(dx * dx + dy * dy);
        }
    }
    return 2.0 / static_cast<double>(N * (N - 1)) * distance;
}

double avg_pair_dist(const CanonicalEnsemble<Particle3D>::State &state) {
    const auto N = state.size();
    double distance = 0.0;
    for (auto it = state.cbegin(), end = state.cend(); it != end; ++it) {
        for (auto it2 = state.cbegin(); it2 != it; ++it2) {
            const auto dx = it->x - it2->x;
            const auto dy = it->y - it2->y;
            const auto dz = it->z - it2->z;
            distance += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }
    return 2.0 / static_cast<double>(N * (N - 1)) * distance;
}
