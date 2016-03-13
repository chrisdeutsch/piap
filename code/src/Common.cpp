#include "Common.h"

std::function<Particle2D(const Particle2D &, std::mt19937 &)>
make_proposal_function(double stddev, double box_length) {
    auto limit = box_length / 2.0;
    std::normal_distribution<> dist(0.0, stddev);

    auto proposal_function = [dist, limit](const Particle2D &p,
                                           std::mt19937 &rng) mutable {
        auto ret = p;

        do {
            ret.x = p.x + dist(rng);
        } while (ret.x < -limit || ret.x > limit);

        do {
            ret.y = p.y + dist(rng);
        } while (ret.y < -limit || ret.y > limit);

        return ret;
    };

    return proposal_function;
}

std::function<Particle2D(const Particle2D &, std::mt19937 &)>
unif_proposal_function(double delta, double box_length) {
    auto limit = box_length / 2.0;
    std::uniform_real_distribution<> unif_dist(-delta, delta);

    auto proposal_function = [unif_dist, limit](const Particle2D &p,
                                         std::mt19937 &rng) mutable {
        auto ret = p;

        do {
            ret.x = p.x + unif_dist(rng);
        } while (ret.x < -limit || ret.x > limit);

        do {
            ret.y = p.y + unif_dist(rng);
        } while (ret.y < -limit || ret.y > limit);

        return ret;
    };

    return proposal_function;
}

std::vector<Particle2D> random_state(double box_length, unsigned pair_num) {
    thread_local std::mt19937 rng(std::random_device{}());

    std::vector<Particle2D> ret;
    std::uniform_real_distribution<> unif(-box_length / 2.0, box_length / 2.0);

    for (unsigned i = 0; i < pair_num; ++i) {
        ret.emplace_back(1.0, unif(rng), unif(rng));
        ret.emplace_back(-1.0, unif(rng), unif(rng));
    }
    return ret;
}
