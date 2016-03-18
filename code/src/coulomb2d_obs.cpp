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
    const auto init_state = random_state(width, n_pairs);
    auto prop_func = unif_proposal_function(sigma, width);

    using Ensemble = CanonicalEnsemble<Particle2D>;
    using PotentialPtr = double (*)(const Particle2D &, const Particle2D &);

    Ensemble ensemble(init_state, beta, static_cast<PotentialPtr>(coulomb_core),
                      prop_func);

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
        if (beta >= 1.0 && beta < 5.0) {
            return 2.0;
        }
        else if (beta >= 5.0 && beta < 10.0) {
            return 2.0 + (1.4 - 2.0) / (10.0 - 5.0) * (beta - 5.0);
        }
        else if (beta >= 10.0 && beta < 20.0) {
            return 1.4 + (0.85 - 1.4) / (20.0 - 10.0) * (beta - 10.0);
        }
        else if (beta >= 20.0 && beta < 40.0) {
            return 0.85 + (0.5 - 0.85) / (40.0 - 20.0) * (beta - 20.0);
        }
        else if (beta >= 40.0 && beta < 80.0) {
            return 0.5 + (0.3 - 0.5) / (80.0 - 40.0) * (beta - 40.0);
        }
        else if (beta >= 80.0 && beta < 160.0) {
            return 0.3 + (0.22 - 0.3) / (160.0 - 80.0) * (beta - 80.0);
        }
        else if (beta >= 160.0 && beta < 220.0) {
            return 0.22 + (0.18 - 0.22) / (220.0 - 160.0) * (beta - 160.0);
        }
        else if (beta >= 220.0 && beta < 380.0) {
            return 0.18 + (0.14 - 0.18) / (380.0 - 220.0) * (beta - 220.0);
        }
        else if (beta >= 380.0 && beta < 500.0) {
            return 0.14 + (0.12 - 0.14) / (500.0 - 380.0) * (beta - 380.0);
        }
        else {
            return 0.12;
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
        double beta = 1.0;
        while (beta < 500.0) {
            std::cout << "beta: " << beta << std::endl;
            double avg_pair_d, acceptance;

            for (int i = 0; i < 5; ++i) {
                auto calc1 = std::async(std::launch::async, calc_obs, beta,
                    gauge_curve_unif_30(beta), 15.0, 20, 1000000);
                auto calc2 = std::async(std::launch::async, calc_obs, beta,
                    gauge_curve_unif_30(beta), 15.0, 20, 1000000);
                auto calc3 = std::async(std::launch::async, calc_obs, beta,
                    gauge_curve_unif_30(beta), 15.0, 20, 1000000);

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
