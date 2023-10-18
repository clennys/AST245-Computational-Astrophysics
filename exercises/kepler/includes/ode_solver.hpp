#pragma once
#include "particle.hpp"
class ODESolver {
  public:
    enum class ODEScheme {
        RungeKutta2nd,
        RungeKutta4th,
        LeapFrop,
        SemiImplEuler,
        ExplicitEuler
    };
    enum class ODEDerivatives { KeplerianOrbits };
    /// Enum to hold transformation type
    enum class TransElemem { Position, Velocity };
    // TODO: (aver) add enum class for what kind of system

    ODESolver(ODEScheme scheme, ODEDerivatives derivative_function,
              double timesteps);
    ODESolver(ODESolver &&) = default;
    ODESolver(const ODESolver &) = default;
    ODESolver &operator=(ODESolver &&) = delete;
    ODESolver &operator=(const ODESolver &) = delete;
    ~ODESolver();

    /// @brief Solve a an ODE system defined by the constructor
    auto solve_system(Particle &init_particle, const size_t &period)
        -> std::vector<Particle>;
    // ============================================================================================
    // Public Static Helper Functions
    // ============================================================================================
    static auto transform_vec2d(const std::vector<Particle> &particles,
                                const TransElemem &type)
        -> std::pair<std::vector<double>, std::vector<double>>;

  private:
    // ============================================================================================
    // member variables
    // ============================================================================================
    ODEScheme m_scheme;
    ODEDerivatives m_derivFunction;

    // ============================================================================================
    // constants
    // ============================================================================================

    static constexpr const int k_GM = 1;
    // static constexpr const int k_particle_mass = 1;
    const double k_diff_step;
    /// eccentricity 0 <= e < 1
    static constexpr const double k_eccentricity = 0.001;
    /// semimajor axis a
    // static constexpr const double k_semimajor_axis = 1 / (1 -
    // k_eccentricity);

    // ============================================================================================
    // Helper Methods
    // ============================================================================================

    auto derive(const Particle &particle_n, double time_n) -> Particle;
    auto solver_step(const Particle &particle_n, double time_n) -> Particle;

    // ============================================================================================
    // Derivate Methods
    // ============================================================================================

    // TODO: (aver) could just return an new particle
    auto derive_keplerian_orbit(const Particle &particle, double time)
        -> Particle;
    // -> std::pair<Eigen::Vector2d, Eigen::Vector2d>;

    // ============================================================================================
    // Schemes
    // ============================================================================================

    auto explicit_euler(const Particle &particle_n, double time_n) -> Particle;
};
