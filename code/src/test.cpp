#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

#include "RandomWalkMetropolis.h"

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

    auto proposal = [&](const Vec2d &v, Vec2d &dest) {
        dest.x = v.x + dist(gen);
        dest.y = v.y + dist(gen);
    };

    std::vector<Vec2d> samples;
    samples.reserve(1000000);

    Vec2d x0;
    RandomWalkMetropolis<Vec2d> rwmetro(target, proposal, x0);

    bool accepted;
    unsigned int accepted_cnt = 0;
    for (int i = 0; i < 1000000; ++i) {
        std::tie(x0, accepted) = rwmetro.step();
        samples.push_back(x0);
        if (accepted) {
            ++accepted_cnt;
        }
    }

    double sum = 0.0;
    std::ofstream os("out.tsv");
    for (auto it = samples.cbegin(); os && it != samples.cend(); ++it) {
        sum += it->x * it->x + it->y * it->y;
        //os << it->x << "\t" << it->y << "\n";
    }

    std::cout << "<r^2> = " << sum / 1000000.0 << std::endl;
	std::cout << "numerical value = " << 1.0040565 << std::endl;
    std::cout << "acceptance rate = "
              << static_cast<double>(accepted_cnt) / 1000000.0 << std::endl;

    return 0;
}
