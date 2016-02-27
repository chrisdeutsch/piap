#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

/**
 * Random Walk Metropolis
 */
template <typename State>
std::pair<State, bool>
metropolis_update(std::function<double(const State &)> target,
                  std::function<State(const State &)> proposal,
                  const State &state) {
    // 0. Initialize RNG
    static std::mt19937 gen(std::random_device{}());
    static std::uniform_real_distribution<> dist;

    // 1. Propose a new state
    auto proposed = proposal(state);

    // 2. Accept-Reject step
    auto accept_prob = target(proposed) / target(state);

    // TODO: Benchmark multiple branching vs naive method
    auto unif = dist(gen);
    if (unif < accept_prob) {
        return std::make_pair(proposed, true);
    } else {
        return std::make_pair(state, false);
    }
}

/**
 * For testing purposes
 */
struct Vec2d {
    double x = 0.0;
    double y = 0.0;
};

int main() {
    // Some target distribution for testing
    auto target = [](const Vec2d &v) {
        return std::exp(-5.0 * std::abs(v.x * v.x + v.y * v.y - 1));
    };

    std::mt19937 gen(std::random_device{}());
    std::normal_distribution<> dist(0.0, 0.5);

    auto proposal = [&](const Vec2d &v) {
        Vec2d ret = v;
        ret.x += dist(gen);
        ret.y += dist(gen);
        return ret;
    };

    std::vector<Vec2d> samples;
    samples.reserve(100000);

    Vec2d x0;
    bool accepted;
    unsigned int accepted_cnt = 0;
    for (int i = 0; i < 100000; ++i) {
        std::tie(x0, accepted) = metropolis_update<Vec2d>(target, proposal, x0);
        samples.push_back(x0);
        if (accepted) {
            ++accepted_cnt;
        }
    }

    double sum = 0.0;
    std::ofstream os("out.tsv");
    for (auto it = samples.cbegin(); os && it != samples.cend(); ++it) {
        sum += it->x * it->x + it->y * it->y;
        os << it->x << "\t" << it->y << "\n";
    }

    std::cout << "<r^2> = " << sum / 100000.0 << std::endl;
	std::cout << "numerical value = " << 1.0040565 << std::endl;
    std::cout << "acceptance rate = "
              << static_cast<double>(accepted_cnt) / 100000.0 << std::endl;

    return 0;
}
