#ifndef CANONICALENSEMBLE_H_
#define CANONICALENSEMBLE_H_

#include <cmath>
#include <functional>
#include <random>
#include <vector>

/**
 * CanonicalEnsemble
 *
 * Simulates a canonical ensemble at a specified temperature and potential using
 * a Random-Walk Metropolis-Algorithm using a single particle update with a
 * specified proposal distribution.
 */
template <typename ParticleState>
class CanonicalEnsemble {
public:
    /**
     * Type specifying the state of the system.
     */
    using State = std::vector<ParticleState>;

    /**
     * Type specifying the interparticle potential.
     */
    using PotentialFunction =
        std::function<double(const ParticleState &, const ParticleState &)>;

    /**
     * Type of a function proposing a new state for a single particle from the
     * current state.
     */
    using ProposalFunction =
        std::function<ParticleState(const ParticleState &, std::mt19937 &)>;

public:
    /**
     * Constructor taking the initial state of the simulation, thermodynamic
     * beta, interparticle potential and proposal function.
     */
    CanonicalEnsemble(const State &initial_state, double beta,
                      PotentialFunction potential_func,
                      ProposalFunction proposal_func)
        : rng(std::random_device{}()), unif_index(0, initial_state.size() - 1),
          potential_func(potential_func), proposal_func(proposal_func),
          state(initial_state), beta(beta) {}

    /**
     * Returns a reference to the current state of the simulation.
     */
    const State &get_state() const { return state; }

    /**
     * Executes a single step of the Random-Walk Metropolis-Algorithm
     *
     * returns: bool - indicating whether the step was accepted.
     */
    bool step() {
        // 1. Choose particle at random and calculate its potential
        const auto idx = unif_index(rng);
        const auto current_pot = potential(idx);

        // 2. Propose a new state for the chosen particle and calculate the new
        //    potential
        const auto current = state[idx];
        state[idx] = proposal_func(current, rng);
        const auto proposed_pot = potential(idx);

        // 3. Accept-reject step
        const auto accept_prob = std::exp(-beta * (proposed_pot - current_pot));

        const auto accepted = unif_real(rng) < accept_prob;
        if (!accepted) {
            // Revert proposal
            state[idx] = current;
        }
        return accepted;
    }

private:
    /**
     * Calculates the potential of the specified particle
     */
    double potential(typename State::size_type ix) const {
        const auto p = state.cbegin() + ix;
        auto potential_energy = 0.0;

        for (auto it = state.cbegin(), end = state.cend(); it != end; ++it) {
            if (it != p) {
                potential_energy += potential_func(*it, *p);
            }
        }
        return potential_energy;
    }

private:
    std::mt19937 rng;
    std::uniform_real_distribution<double> unif_real;
    std::uniform_int_distribution<typename State::size_type> unif_index;

    PotentialFunction potential_func;
    ProposalFunction proposal_func;

    State state;

    double beta;
};

#endif // CANONICALENSEMBLE_H_
