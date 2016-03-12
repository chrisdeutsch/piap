#include <cstddef>
#include <fstream>
#include <future>
#include <iostream>
#include <string>

#include "Common.h"

/**
 * Calculates average pair distance from samples
 * returns: pair containing expectation value and acceptance rate
 */
std::pair<double, double> calc_obs(double beta, double sigma, double width,
                                   std::size_t n_pairs, std::size_t n_samples) {
    const auto init_state = random_state(width, n_pairs);
    auto prop_func = make_proposal_function(sigma, width);

    using Ensemble = CanonicalEnsemble<Particle2D>;
    using PotentialPtr = double (*)(const Particle2D &, const Particle2D &);

    Ensemble ensemble(init_state, beta, static_cast<PotentialPtr>(coulomb_core),
                      prop_func);

    std::size_t acceptance_cnt = 0;
    double expect_val = 0.0;

    // Burn-in (TODO: make variable)
    for (std::size_t i = 0; i < 1000; ++i) {
        ensemble.step();
    }

    // Calculation of the expectation value
    for (std::size_t i = 0; i < n_samples; ++i) {
        if (ensemble.step()) {
            ++acceptance_cnt;
        }
        expect_val += avg_pair_dist(ensemble.get_state());
    }
    return std::make_pair(expect_val / n_samples,
                          static_cast<double>(acceptance_cnt) / n_samples);
}

int main() {
    // First test generating observables
    auto gauge_curve = [](double beta) {
        if (beta >= 1.0 && beta < 10.0) {
            return 3.0 + (1.0 - 3.0) / (10. - 1.) * (beta - 1.0);
        } else if (beta >= 10.0 && beta < 30.0) {
            return 1.0 + (0.4 - 1.0) / (30. - 10.) * (beta - 10.0);
        } else if (beta >= 30.0 && beta < 85.0) {
            return 0.4 + (0.25 - 0.4) / (85. - 30.) * (beta - 30.0);
        } else if (beta >= 85.0 && beta < 100.0) {
            return 0.25 + (0.2 - 0.25) / (100. - 85.) * (beta - 85.0);
        } else if (beta >= 100.0 && beta < 200.0) {
            return 0.2 + (0.15 - 0.2) / (200. - 100.) * (beta - 100.0);
        } else if (beta >= 200.0 && beta < 300.0) {
            return 0.15 + (0.1 - 0.15) / (300. - 200.) * (beta - 200.0);
        } else {
            return 0.1;
        }
    };

    std::ofstream os("observables.csv");

    os << "beta,acc,obs\n";

    double beta = 1.0;
    while (beta < 400.0) {
        std::cout << "beta: " << beta << std::endl;
        double avg_pair_d, acceptance;

        for (int i = 0; i < 33; ++i) {
            auto calc1 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve(beta), 15.0, 20, 10000);
            auto calc2 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve(beta), 15.0, 20, 10000);
            auto calc3 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve(beta), 15.0, 20, 10000);

            std::tie(avg_pair_d, acceptance) = calc1.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";

            std::tie(avg_pair_d, acceptance) = calc2.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";

            std::tie(avg_pair_d, acceptance) = calc3.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";
        }
        beta *= 1.06;
    }

    return 0;
}
