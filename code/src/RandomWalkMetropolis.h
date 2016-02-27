#include <functional>
#include <random>
#include <utility>

/**
 * RandomWalkMetropolis
 *
 * Implementation of a generic Random-Walk Metropolis-Algorithm with customized
 * target distribution and proposal function.
 */
template <typename State>
class RandomWalkMetropolis {
public:
    /**
     * Function proportional to the desired probability distribution
     */
    using TargetDistribution = std::function<double(const State &)>;

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
     * Constructor taking the target distribution, proposal function and an
     * initial state.
     */
    RandomWalkMetropolis(TargetDistribution target_distribution,
                         ProposalFunction proposal_function,
                         const State &initial_state)
        : rng(std::random_device{}()), target_distribution(target_distribution),
          proposal_function(proposal_function), current_state(initial_state),
          proposed_state(initial_state),
          weight_current_state(target_distribution(initial_state)) {}

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
        const auto weight_proposed_state = target_distribution(proposed_state);
        const auto accept_prob = weight_proposed_state / weight_current_state;

        // unif in interval [0.0, 1.0)
        const auto unif = unif_dist(rng);
        if (unif < accept_prob) {
            // accept proposed state
            std::swap(current_state, proposed_state);
            weight_current_state = weight_proposed_state;
            return std::make_pair(current_state, true);
        } else {
            // reject proposed state
            return std::make_pair(current_state, false);
        }
    }

private:
    std::mt19937 rng;
    std::uniform_real_distribution<> unif_dist;

    TargetDistribution target_distribution;
    ProposalFunction proposal_function;

    State current_state;
    State proposed_state;

    /**
     * Since the evaluation of the target distribution can be expensive, the
     * weight of the current state is cached.
     */
    double weight_current_state;
};
