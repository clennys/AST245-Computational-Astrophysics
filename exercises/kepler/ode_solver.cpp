#include "ode_solver.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <numbers>
#include <utility>
#include <vector>

#include "particle.hpp"

ODESolver::ODESolver(ODEScheme scheme, ODEDerivatives derivative_function,
                     double timesteps)
    : m_scheme(scheme), m_derivFunction(derivative_function),
      k_diff_step(timesteps) {}

auto ODESolver::solve_system(Particle &init_particle, const size_t &period)
    -> std::vector<Particle> {

    const auto k_orbital_period =
        (2 * std::numbers::pi) / std::pow(1 - k_eccentricity, 1.5);
    const auto k_final_time = k_orbital_period * static_cast<double>(period);
    std::cout << "k_final_time: " << k_final_time << "\n";
    std::cout << "k_diff_step: " << k_diff_step << "\n";
    auto current_time = 0.;
    std::vector<Particle> particles;

    // setup particle
    init_particle.compute_energy();
    init_particle.compute_eccentricity();
    init_particle.compute_angular_momentum();

    while (current_time < k_final_time) {

        auto new_particle = solver_step(init_particle, current_time);
        // std::cerr << "DEBUGPRINT[1]: kepler.cpp:32: particles="
        //           << particles.size() << std::endl;
        new_particle.compute_angular_momentum();
        new_particle.compute_eccentricity();
        new_particle.compute_energy();

        particles.push_back(new_particle);

        current_time += this->k_diff_step;
        // std::cout << "current_time:" << current_time << "\n";
    }

    return particles;
}

ODESolver::~ODESolver() {}

// ============================================================================================
// Derivate Method Implementations
// ============================================================================================

auto ODESolver::derive(const Particle &particle_n, double time_n) -> Particle {
    switch (m_derivFunction) {
    case ODEDerivatives::KeplerianOrbits:
    default:
        return derive_keplerian_orbit(particle_n, time_n);
    }
}
auto ODESolver::derive_keplerian_orbit(const Particle &particle, double time)
    -> Particle {
    const auto new_position =
        (-k_GM / std::pow(particle.position.norm(), 3)) * particle.position;

    return Particle(new_position, particle.velocity, particle.k_mass);
}

// ============================================================================================
// Scheme Implementations
// ============================================================================================

auto ODESolver::solver_step(const Particle &particle_n, double time_n)
    -> Particle {
    switch (m_scheme) {
    case ODEScheme::RungeKutta2nd:
    case ODEScheme::RungeKutta4th:
    case ODEScheme::LeapFrop:
    case ODEScheme::SemiImplEuler:
    case ODEScheme::ExplicitEuler:
        return explicit_euler(particle_n, time_n);
    default:
        break;
    }
}
auto ODESolver::explicit_euler(const Particle &particle_n, double time_n)
    -> Particle {
    auto der_particle = derive(particle_n, time_n);
    return Particle(
        der_particle.position * this->k_diff_step + particle_n.position,
        der_particle.velocity * this->k_diff_step + particle_n.velocity,
        der_particle.k_mass);
}
