#ifndef CANONICALENSEMBLE_H_
#define CANONICALENSEMBLE_H_

#include <cmath>
#include <functional>
#include <random>
#include <utility>

/**
 * CanonicalEnsemble
 *
 * Simulation of a canonical ensemble at a specified temperature using a
 * Random-Walk Metropolis-Algorithm with customized Hamiltonian and proposal
 * function.
 */
template <typename State>
class CanonicalEnsemble {
public:
    /**
     * Hamiltonian of the system as a function of the realized state.
     */
    using Hamiltonian = std::function<double(const State &)>;

    /**
     * Function proposing a new state from the current state "current"
     * and saves it into "destination". "destination" is guaranteed to be a in a
     * valid state.
     *
     * Signature: void(const State &current, State &destination)
     */
    using ProposalFunction = std::function<void(const State &, State &)>;

public:
    /**
     * Constructor taking the system's Hamiltonian, the thermodynamic beta,
     * a proposal function and the initial state of the simulation.
     */
    CanonicalEnsemble(Hamiltonian hamiltonian, double beta,
                      ProposalFunction proposal_function,
                      const State &initial_state)
        : rng(std::random_device{}()), hamiltonian(hamiltonian),
          proposal_function(proposal_function), current_state(initial_state),
          proposed_state(initial_state), beta(beta),
          current_energy(hamiltonian(initial_state)) {}

    /**
     * Evaluates a single step of the Random-Walk Metropolis-Algorithm
     *
     * Returns a pair of the sampled state and a boolean indicating whether a
     * new state was accepted.
     */
    std::pair<State, bool> step() {
        // 1. Propose a new state
        proposal_function(current_state, proposed_state);

        // 2. Accept-reject step
        const auto proposed_energy = hamiltonian(proposed_state);
        const auto accept_prob =
            std::exp(-beta * (proposed_energy - current_energy));

        // unif_dist: uniformly distributed in [0.0, 1.0)
        const auto accepted = unif_dist(rng) < accept_prob;
        if (accepted) {
            std::swap(current_state, proposed_state);
            current_energy = proposed_energy;
        }
        return std::make_pair(current_state, accepted);
    }

private:
    std::mt19937 rng;
    std::uniform_real_distribution<> unif_dist;

    Hamiltonian hamiltonian;
    ProposalFunction proposal_function;

    State current_state;
    State proposed_state;

    double beta;

    /**
     * Since the evaluation of the Hamiltonian can be expensive, the energy of
     * the current state is cached.
     */
    double current_energy;
};

#endif // CANONICALENSEMBLE_H_
