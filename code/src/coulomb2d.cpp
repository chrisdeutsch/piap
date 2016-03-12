#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>

#include "Common.h"

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

    // Convenience typedefs
    using Ensemble = CanonicalEnsemble<Particle2D>;
    using PotentialPtr = double (*)(const Particle2D &, const Particle2D &);

    auto initial_state = random_state(side_length, particle_num);
    auto proposal_function =
        make_proposal_function(proposal_stddev, side_length);

    // Simulate samples of the canonical ensemble using the Random-Walk
    // Metropolis-Algorithm
    Ensemble ensemble(initial_state, beta,
                      static_cast<PotentialPtr>(coulomb_core),
                      proposal_function);

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
