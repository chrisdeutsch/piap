#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "CanonicalEnsemble.h"

struct Particle2D {
    Particle2D(double q, double x, double y) : q(q), x(x), y(y) {}
    double q = 1.0;
    double x = 0.0;
    double y = 0.0;
};

using ParticleState = std::vector<Particle2D>;

double
hamiltonian_w_pot(const ParticleState &state,
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
    // Standard deviation of truncated normal distribution in proposal function
    double proposal_stddev = 0.01;
    // Side length of the box
    double side_length = 12.0;
    // Number of oppositely charged particle pairs
    unsigned particle_num = 25;
    // Thermodynamic beta
    double beta = 200.0;
    // Number of samples
    std::size_t sample_num = 10000;
    // Output
    bool save_output = true;
    std::string filename = "out.tsv";

    // Interparticle potential
    auto coulomb_w_core = [](const Particle2D &a, const Particle2D &b) {
        const auto dx = a.x - b.x;
        const auto dy = a.y - b.y;
        const auto distance = std::sqrt(dx * dx + dy * dy);
        return a.q * b.q / distance + std::pow(distance, -8.0);
    };

    // Hamiltonian with Coulomb potential and hard cores
    using namespace std::placeholders;
    auto hamiltonian = std::bind(hamiltonian_w_pot, _1, coulomb_w_core);

    // Proposal function
    auto proposal_fun = [=](const ParticleState &current,
                            ParticleState &destination) {
        auto curr_it = current.cbegin(), curr_end = current.cend();
        auto dest_it = destination.begin();

        while (curr_it != curr_end) {
            dest_it->q = curr_it->q;
            dest_it->x = trunc_normal(curr_it->x, proposal_stddev,
                                      -side_length / 2.0, side_length / 2.0);
            dest_it->y = trunc_normal(curr_it->y, proposal_stddev,
                                      -side_length / 2.0, side_length / 2.0);

            ++curr_it;
            ++dest_it;
        }
    };

    // Randomize initial state
    ParticleState initial_state;

    std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<> unif(-side_length / 2.0,
                                          side_length / 2.0);

    for (unsigned n = 0; n < particle_num; ++n) {
        initial_state.emplace_back(1.0, unif(rng), unif(rng));
        initial_state.emplace_back(-1.0, unif(rng), unif(rng));
    }

    // Simulate samples of the canonical ensemble using the Random-Walk
    // Metropolis-Algorithm
    CanonicalEnsemble<ParticleState> rw_metro(hamiltonian, beta, proposal_fun,
                                              initial_state);
    std::vector<ParticleState> samples;
    samples.reserve(sample_num);

    ParticleState state;
    bool accepted;
    std::size_t accepted_cnt = 0;

    for (std::size_t i = 0; i < sample_num; ++i) {
        std::tie(state, accepted) = rw_metro.step();
        samples.push_back(std::move(state));

        if (accepted) {
            ++accepted_cnt;
        }
    }

    std::cout << "Acceptance probability: "
              << static_cast<double>(accepted_cnt) / sample_num << std::endl;

    // Save output to .tsv
    if (save_output) {
        std::ofstream os(filename);
        for (auto it = samples.cbegin(), end = samples.cend(); os && it != end;
             ++it) {
            for (const auto &particle : *it) {
                os << particle.q << "\t" << particle.x << "\t" << particle.y
                   << "\t";
            }
            os << "\n";
        }
    }

    return 0;
}
