#include "particle.hpp"
#include "logging.hpp"

#include <algorithm>
#include <execution>
#include <format>
#include <iostream>
#include <mutex>

#include <omp.h>

double Particles::g_total_mass = 0.;

auto Particle3D::calc_orign_distance() -> void { this->distance = this->position.norm(); }

auto Particle3D::print_summary() const -> void {
    std::cout << "Position = " << this->position << std::endl;
    std::cout << "Potential = " << this->potential << std::endl;
    std::cout << "Velocity = " << this->velocity << std::endl;
    std::cout << "Distance = " << this->distance << std::endl;
}

auto Particles::get_max_distance(const std::vector<Particle3D> &particles) -> Particle3D {
    auto result =
        std::max_element(particles.begin(), particles.end(), [](Particle3D a, Particle3D b) {
            return a.distance < b.distance;
        });

    result->print_summary();

    return *result;
}
auto Particles::calc_total_mass(const std::vector<Particle3D> &particles) -> void {
    g_total_mass = 0.;

#pragma omp parallel for
    for (const auto &part : particles) {
#pragma omp atomic
        g_total_mass += part.mass;
    }

    // NOTE: (aver) consider simpler omp loop
    // std::mutex mtx;
    // std::for_each(
    //     std::execution::par, particles.begin(), particles.end(), [&mtx](const Particle3D &part) {
    //         std::lock_guard<std::mutex> guard(mtx); // Protect m_shells
    //         Particles::g_total_mass += part.mass;
    //     });

    Logging::info(std::format("Total mass of system: {}", Particles::g_total_mass));
}
