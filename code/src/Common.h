#ifndef COMMON_H_
#define COMMON_H_

#include <functional>
#include "Particle.h"

/**
 * Closure to create a proposal function with the given parameters.
 */
std::function<Particle2D(const Particle2D &, std::mt19937 &)>
make_proposal_function(double stddev, double box_length);

/**
 * Uniform out of box with side length 2*delta
 */
std::function<Particle2D(const Particle2D &, std::mt19937 &)>
unif_proposal_function(double delta, double box_length);

/**
 * Creates a state with uniformly distributed particles.
 */
std::vector<Particle2D> random_state(double box_length, unsigned pair_num);

#endif // COMMON_H_
