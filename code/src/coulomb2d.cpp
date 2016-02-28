#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "RandomWalkMetropolis.h"

struct Particle2D {
    Particle2D(double q, double x, double y) : q(q), x(x), y(y) {}
    double q = 1.0;
    double x = 0.0;
    double y = 0.0;
};

using ParticleState = std::vector<Particle2D>;

double tot_pot_energy(const ParticleState &state,
                      std::function<double(const Particle2D &a,
                                           const Particle2D &b)> pot_energy) {
    double epot = 0.0;
    for (auto it = state.cbegin(), end = state.cend(); it != end; ++it) {
        for (auto it2 = state.cbegin(); it2 != it; ++it2) {
            epot += pot_energy(*it, *it2);
        }
    }
    return epot;
}

// Put this in a class?
double trunc_normal(double mean, double stddev, double lower, double upper) {
    static std::mt19937 rng(std::random_device{}());
    std::normal_distribution<> dist(mean, stddev);

    auto rand = dist(rng);
    while (rand < lower || rand > upper) {
        rand = dist(rng);
    }

    return rand;
}

int main() {
    // Setup initial state
    ParticleState initial_state;

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> unif(-10.0, 10.0);
    for (unsigned n = 0; n < 20; ++n) {
        initial_state.emplace_back(1.0, unif(rng), unif(rng));
        initial_state.emplace_back(-1.0, unif(rng), unif(rng));
    }

    // Interparticle potential
    auto coulomb_w_core = [](const Particle2D &a, const Particle2D &b) {
        const auto dx = a.x - b.x;
        const auto dy = a.y - b.y;
        const auto distance = std::sqrt(dx * dx + dy * dy);
        return a.q * b.q / distance + std::pow(distance, -8.0);
    };

    // Target distribution
    // FIX: this distribution leads to infinities (RandomWalkMetropolis ->
    //      CanonicalEnsemble)
    //      use delta_E instead of exp(-E_proposed / T) / exp(-E_current / T)
    //      the large exponents in the exponential lead to infinities, since
    //      exp(-large_val) approx 0
    double temp = 1.0;
    auto target_dist = [=](const ParticleState &state) {
        return std::exp(-tot_pot_energy(state, coulomb_w_core) / temp);
    };

    // Proposal function
    auto proposal_fun = [&](const ParticleState &current,
                            ParticleState &destination) {
        auto curr_it = current.cbegin(), curr_end = current.cend();
        auto dest_it = destination.begin(), dest_end = destination.end();

        while (curr_it != curr_end) {
            dest_it->q = curr_it->q;
            dest_it->x = trunc_normal(curr_it->x, 0.2, -10.0, 10.0);
            dest_it->y = trunc_normal(curr_it->y, 0.2, -10.0, 10.0);

            ++curr_it;
            ++dest_it;
        }
    };

    // Random-Walk Metropolis
    RandomWalkMetropolis<ParticleState> rw_metro(target_dist, proposal_fun,
                                                 initial_state);

    std::vector<ParticleState> samples;
    samples.reserve(10000);

    int acceptance = 0;

    ParticleState state;
    bool accepted;
    for (int i = 0; i < 10000; ++i) {
        std::tie(state, accepted) = rw_metro.step();
        samples.push_back(std::move(state));

        if (accepted) {
            ++acceptance;
        }
    }

    std::cout << "acceptance rate: "
              << static_cast<double>(acceptance) / 10000.0 << std::endl;
    std::cout << "accepted: " << acceptance << std::endl;

    std::ofstream os("out.tsv");
    auto it = samples.cbegin(), end = samples.cend();
    while (os && it != end) {
        for (const auto &particle : *it) {
            os << particle.q << "\t" << particle.x << "\t" << particle.y
               << "\t";
        }
        os << "\n";
        ++it;
    }

    return 0;
}
