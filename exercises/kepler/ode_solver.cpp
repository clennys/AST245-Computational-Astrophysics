#include "ode_solver.hpp"

#include <Eigen/Dense>
#include <future>
#include <iostream>
#include <numbers>
#include <ranges>
#include <utility>
#include <vector>

#include "particle.hpp"

auto ODESolver::transform_vec2d(const std::vector<Particle> &particles,
                                const TransElemem &type)
    -> std::pair<std::vector<double>, std::vector<double>> const {

    auto future1 = std::async(std::launch::async, [&]() {
        auto transformed_x =
            particles |
            std::views::transform([&type](const Particle &particle) {
                if (type == TransElemem::Position)
                    return particle.position.x();
                if (type == TransElemem::Velocity)
                    return particle.velocity.x();
            });
        std::vector<double> transformed_vector_x(transformed_x.begin(),
                                                 transformed_x.end());
        return transformed_vector_x;
    });
    auto future2 = std::async(std::launch::async, [&]() {
        auto transformed_y =
            particles |
            std::views::transform([&type](const Particle &particle) {
                if (type == TransElemem::Position)
                    return particle.position.y();
                if (type == TransElemem::Velocity)
                    return particle.velocity.y();
            });
        std::vector<double> transformed_vector_y(transformed_y.begin(),
                                                 transformed_y.end());
        return transformed_vector_y;
    });
    // we need to transform the vectors of particles into positions for
    // plotting
    return std::pair(future1.get(), future2.get());
}

ODESolver::ODESolver(ODEScheme scheme, ODEDerivatives derivative_function,
                     double timesteps)
    : m_scheme(scheme), m_derivFunction(derivative_function),
      diff_step(timesteps) {}

auto ODESolver::solve_system(Particle &init_particle, const size_t &period,
                             const double &k_eccentricity)
    -> std::vector<Particle> {

    const auto k_orbital_period =
        (2 * std::numbers::pi) / std::pow(1 - k_eccentricity, 1.5);
    const auto k_final_time = k_orbital_period * static_cast<double>(period);
    auto current_time = 0.;
    std::vector<Particle> particles;

    // setup particle
    init_particle.compute_energy();
    init_particle.compute_eccentricity();
    init_particle.compute_angular_momentum();

    while (current_time < k_final_time) {
        auto new_particle = solver_step(init_particle, current_time);
        new_particle.compute_angular_momentum();
        new_particle.compute_eccentricity();
        new_particle.compute_energy();

        particles.push_back(new_particle);

        // update particle and time
        init_particle = new_particle;
        current_time += this->diff_step;
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
    (void)time;
    const auto new_position =
        (-this->k_GM / std::pow(particle.position.norm(), 3)) *
        particle.position;
    return Particle(particle.velocity, new_position, particle.k_mass);
}

// ============================================================================================
// Scheme Implementations
// ============================================================================================

auto ODESolver::solver_step(const Particle &particle_n, double time_n)
    -> Particle {
    switch (m_scheme) {
    case ODEScheme::RungeKutta2nd:
        return runge_kutta_2nd_order(particle_n, time_n);
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
        der_particle.position * this->diff_step + particle_n.position,
        der_particle.velocity * this->diff_step + particle_n.velocity,
        der_particle.k_mass);
}
auto ODESolver::runge_kutta_2nd_order(const Particle &particle_n, double time_n)
    -> Particle {
    auto particle_1 = derive(particle_n, time_n);
    auto particle_2 = derive(
        Particle(particle_n.position + particle_1.position * this->diff_step,
                 particle_n.velocity + particle_1.velocity * this->diff_step,
                 particle_n.k_mass),
        time_n);
    auto derived_particle = Particle(
        (particle_1.position + particle_2.position) / 2,
        (particle_1.velocity + particle_2.velocity) / 2, particle_n.k_mass);

    return Particle(
        particle_n.position + derived_particle.position * this->diff_step,
        particle_n.velocity + derived_particle.velocity * this->diff_step,
        particle_n.k_mass);
}
