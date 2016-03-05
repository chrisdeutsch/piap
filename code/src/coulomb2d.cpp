#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include "CanonicalEnsemble.h"

struct Particle2D {
    Particle2D(double q, double x, double y) : q(q), x(x), y(y) {}
    double q = 1.0;
    double x = 0.0;
    double y = 0.0;
};

int main() {
    // Standard deviation of truncated normal distribution in proposal function
    double proposal_stddev = 0.1;
    // Side length of the box
    double side_length = 15.0;
    // Number of oppositely charged particle pairs
    unsigned particle_num = 20;
    // Thermodynamic beta
    double beta = 300.0;
    // Number of samples
    std::size_t sample_num = 100000;
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

    // Proposal Function
    auto proposal_function = [side_length, proposal_stddev](
        const Particle2D &particle, std::mt19937 &rng) {
        auto ret = particle;
        const auto upper = side_length / 2.0;
        const auto lower = -upper;

        std::normal_distribution<> dist_x(particle.x, proposal_stddev);
        std::normal_distribution<> dist_y(particle.y, proposal_stddev);

        do {
            ret.x = dist_x(rng);
        } while (ret.x < lower || ret.x > upper);

        do {
            ret.y = dist_y(rng);
        } while (ret.y < lower || ret.y > upper);

        return ret;
    };

    // Randomize initial state
    using Ensemble = CanonicalEnsemble<Particle2D>;

    std::mt19937 rng(std::random_device{}());
    Ensemble::State initial_state;
    std::uniform_real_distribution<> unif(-side_length / 2.0,
                                          side_length / 2.0);

    for (unsigned n = 0; n < particle_num; ++n) {
        initial_state.emplace_back(1.0, unif(rng), unif(rng));
        initial_state.emplace_back(-1.0, unif(rng), unif(rng));
    }

    // Simulate samples of the canonical ensemble using the Random-Walk
    // Metropolis-Algorithm
    Ensemble ensemble(initial_state, beta, coulomb_w_core, proposal_function);
    std::vector<Ensemble::State> samples;
    samples.reserve(sample_num);

    std::size_t accepted_cnt = 0;
    for (std::size_t i = 0; i < sample_num; ++i) {
        if (ensemble.step()) {
            ++accepted_cnt;
        }
        samples.push_back(ensemble.get_state());
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
