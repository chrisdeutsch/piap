#include <chrono>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "Common.h"

/**
 * Calculates average pair distance from samples
 * returns: pair containing expectation value and acceptance rate
 */
std::pair<double, double> calc_obs(double beta, double sigma, double width,
                                   std::size_t n_pairs, std::size_t n_samples) {
    const auto init_state = random_state_3d(width, n_pairs);
    auto prop_func = unif_proposal_function_3d(sigma, width);

    using Ensemble = CanonicalEnsemble<Particle3D>;
    using PotentialPtr = double (*)(const Particle3D &, const Particle3D &);

    Ensemble ensemble(init_state, beta,
                      static_cast<PotentialPtr>(lennard_jones), prop_func);

    std::size_t acceptance_cnt = 0;
    std::size_t expect_cnt = 0;
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
        if (i % 20 == 0) {
            ++expect_cnt;
            expect_val += avg_pair_dist(ensemble.get_state());
        }
    }
    return std::make_pair(expect_val / expect_cnt,
                          static_cast<double>(acceptance_cnt) / n_samples);
}

int main() {
    auto gauge_curve_unif_30 = [](double beta) {
        if (beta >= 0.1 && beta < 1.0) {
            return 2.0 + (1.8 - 2.0) / (1.0 - 0.1) * (beta - 0.1);
        }
        else if (beta >= 1.0 && beta < 5.0) {
            return 1.8 + (0.3 - 1.8) / (5.0 - 1.0) * (beta - 1.0);
        }
        else if (beta >= 5.0 && beta < 10.0) {
            return 0.3 + (0.15 - 0.3) / (10.0 - 5.0) * (beta - 5.0);
        }
        else if (beta >= 10.0 && beta < 20.0) {
            return 0.15 + (0.09 - 0.15) / (20.0 - 10.0) * (beta - 10.0);
        }
        else if (beta >= 20.0 && beta < 40.0) {
            return 0.09 + (0.06 - 0.09) / (40.0 - 20.0) * (beta - 20.0);
        }
        else if (beta >= 40.0 && beta < 60.0) {
            return 0.06 + (0.05 - 0.06) / (60.0 - 40.0) * (beta - 40.0);
        }
        else if (beta >= 60.0 && beta < 100.0) {
            return 0.05 + (0.038 - 0.05) / (100.0 - 60.0) * (beta - 60.0);
        }
        else {
            return 0.038;
        }
    };

    while (true) {
        // Timing
        const auto start = std::chrono::high_resolution_clock::now();

        // Get timestamp
        const auto time = std::time(nullptr);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&time), "%Y_%m_%d_%H_%M_%S");
        ss << ".csv";

        std::ofstream os(ss.str());

        // Table header
        os << "# Start of simulation: " << std::ctime(&time);
        os << "beta,acc,obs\n" << std::flush;

        // Simulation
        double beta = 0.1;
        while (beta < 100.0) {
            std::cout << "beta: " << beta << std::endl;
            double avg_pair_d, acceptance;

            for (int i = 0; i < 5; ++i) {
                auto calc1 =
                    std::async(std::launch::async, calc_obs, beta,
                               gauge_curve_unif_30(beta), 5.0, 20, 200000);
                auto calc2 =
                    std::async(std::launch::async, calc_obs, beta,
                               gauge_curve_unif_30(beta), 5.0, 20, 200000);
                auto calc3 =
                    std::async(std::launch::async, calc_obs, beta,
                               gauge_curve_unif_30(beta), 5.0, 20, 200000);

                std::tie(avg_pair_d, acceptance) = calc1.get();
                os << beta << "," << acceptance << "," << avg_pair_d << "\n";

                std::tie(avg_pair_d, acceptance) = calc2.get();
                os << beta << "," << acceptance << "," << avg_pair_d << "\n";

                std::tie(avg_pair_d, acceptance) = calc3.get();
                os << beta << "," << acceptance << "," << avg_pair_d << "\n";

                os << std::flush;
            }
            beta *= 1.04;
        }

        const auto end = std::chrono::high_resolution_clock::now();
        const auto elapsed =
            std::chrono::duration_cast<std::chrono::minutes>(end - start);
        std::cout << "Runtime: " << elapsed.count() << " min" << std::endl;
    }

    return 0;
}
