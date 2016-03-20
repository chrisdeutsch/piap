#ifndef COMMON_H_
#define COMMON_H_

#include <functional>
#include "Particle.h"

/**
 * Uniform proposal function in 2d-box with side length: 2 * delta.
 */
std::function<Particle2D(const Particle2D &, std::mt19937 &)>
unif_proposal_function(double delta, double box_length);

/**
 * Uniform proposal function in 3d-box with side length: 2 * delta.
 */
std::function<Particle3D(const Particle3D &, std::mt19937 &)>
unif_proposal_function_3d(double delta, double box_length);

/**
 * Creates a state with uniformly distributed particles in a 2d-box with side
 * length 2 * delta.
 */
std::vector<Particle2D> random_state(double box_length, unsigned pair_num);

/**
 * Creates a state with uniformly distributed particles in a 3d-box with side
 * length 2 * delta.
 */
std::vector<Particle3D> random_state_3d(double box_length, unsigned pair_num);

#endif // COMMON_H_
