#include <chrono>
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

    auto gauge_curve_unif = [](double beta) {
        if (beta >= 1.0 && beta < 2.0) {
            return 2.0 + (1.65 - 2.0) / (2. - 1.) * (beta - 1.0);
        } else if (beta >= 2.0 && beta < 3.0) {
            return 1.65 + (1.4 - 1.65) / (3. - 2.) * (beta - 2.0);
        } else if (beta >= 3.0 && beta < 5.0) {
            return 1.4 + (1.1 - 1.4) / (5. - 3.) * (beta - 3.0);
        } else if (beta >= 5.0 && beta < 8.0) {
            return 1.1 + (0.85 - 1.1) / (8. - 5.) * (beta - 5.0);
        } else if (beta >= 8.0 && beta < 10.0) {
            return 0.85 + (0.75 - 0.85) / (10. - 8.) * (beta - 8.0);
        } else if (beta >= 10.0 && beta < 15.0) {
            return 0.75 + (0.58 - 0.75) / (15. - 10.) * (beta - 10.0);
        } else if (beta >= 15.0 && beta < 20.0) {
            return 0.58 + (0.47 - 0.58) / (20. - 15.) * (beta - 15.0);
        } else if (beta >= 20.0 && beta < 30.0) {
            return 0.47 + (0.35 - 0.47) / (30. - 20.) * (beta - 20.0);
        } else if (beta >= 30.0 && beta < 40.0) {
            return 0.35 + (0.3 - 0.35) / (40. - 30.) * (beta - 30.0);
        } else if (beta >= 40.0 && beta < 50.0) {
            return 0.3 + (0.26 - 0.3) / (50. - 40.) * (beta - 40.0);
        } else if (beta >= 50.0 && beta < 70.0) {
            return 0.26 + (0.21 - 0.26) / (70. - 50.) * (beta - 50.0);
        } else if (beta >= 70.0 && beta < 90.0) {
            return 0.21 + (0.18 - 0.21) / (90. - 70.) * (beta - 70.0);
        } else if (beta >= 90.0 && beta < 120.0) {
            return 0.18 + (0.15 - 0.18) / (120. - 90.) * (beta - 90.0);
        } else if (beta >= 120.0 && beta < 150.0) {
            return 0.15 + (0.135 - 0.15) / (150. - 120.) * (beta - 120.0);
        } else if (beta >= 150.0 && beta < 180.0) {
            return 0.135 + (0.12 - 0.135) / (180. - 150.) * (beta - 150.0);
        } else if (beta >= 180.0 && beta < 220.0) {
            return 0.12 + (0.11 - 0.12) / (220. - 180.) * (beta - 180.0);
        } else if (beta >= 220.0 && beta < 250.0) {
            return 0.11 + (0.105 - 0.11) / (250. - 220.) * (beta - 220.0);
        } else if (beta >= 250.0 && beta < 300.0) {
            return 0.105 + (0.095 - 0.105) / (300. - 250.) * (beta - 250.0);
        } else if (beta >= 300.0 && beta < 350.0) {
            return 0.095 + (0.09 - 0.095) / (350. - 300.) * (beta - 300.0);
        } else if (beta >= 350.0 && beta < 400.0) {
            return 0.09 + (0.082 - 0.09) / (400. - 350.) * (beta - 350.0);
        } else {
            return 0.082;
        }
    };

    std::ofstream os("observables.csv");

    os << "beta,acc,obs\n";

    const auto start = std::chrono::high_resolution_clock::now();

    double beta = 1.0;
    while (beta < 400.0) {
        std::cout << "beta: " << beta << std::endl;
        double avg_pair_d, acceptance;

        for (int i = 0; i < 20; ++i) {
            auto calc1 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve_unif(beta), 15.0, 20, 300000);
            auto calc2 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve_unif(beta), 15.0, 20, 300000);
            auto calc3 = std::async(std::launch::async, calc_obs, beta,
                                    gauge_curve_unif(beta), 15.0, 20, 300000);

            std::tie(avg_pair_d, acceptance) = calc1.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";

            std::tie(avg_pair_d, acceptance) = calc2.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";

            std::tie(avg_pair_d, acceptance) = calc3.get();
            os << beta << "," << acceptance << "," << avg_pair_d << "\n";
        }
        beta *= 1.06;
    }

    const auto end = std::chrono::high_resolution_clock::now();

    auto elapsed =
        std::chrono::duration_cast<std::chrono::minutes>(end - start);

    std::cout << "Runtime: " << elapsed.count() << " min" << std::endl;

    return 0;
}
