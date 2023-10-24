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
    auto solve_system(Particle &init_particle, const size_t &period,
                      const double &k_eccentricity) -> std::vector<Particle>;
    // ============================================================================================
    // Public Static Helper Functions
    // ============================================================================================
    static auto transform_vec2d(const std::vector<Particle> &particles,
                                const TransElemem &type)
        -> std::pair<std::vector<double>, std::vector<double>> const;

  private:
    // ============================================================================================
    // member variables
    // ============================================================================================
    ODEScheme m_scheme;
    ODEDerivatives m_derivFunction;
    double diff_step;

    // ============================================================================================
    // constants
    // ============================================================================================

    static constexpr const int k_GM = 1;

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

    // ============================================================================================
    // Schemes
    // ============================================================================================

    auto explicit_euler(const Particle &particle_n, double time_n) -> Particle;
    auto runge_kutta_2nd_order(const Particle &particle_n, double time_n)
        -> Particle;
};
