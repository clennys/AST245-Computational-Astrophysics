#include "particle.hpp"
#include "logging.hpp"

#include <iostream>

double Particle3D::s_softening = 0.;

auto Particle3D::calc_origin_dist() -> double { return this->m_position.norm(); }

auto Particle3D::update_origin_dist() -> void { this->m_distance = calc_origin_dist(); }

auto Particle3D::print_summary() const -> void {
    Logging::info("Printing summary:");
    std::cout << "Position =\n" << this->m_position << std::endl;
    std::cout << "Velocity =\n" << this->m_velocity << std::endl;
    std::cout << "Potential = " << this->m_potential << std::endl;
    std::cout << "Distance = " << this->m_distance << std::endl;
    Logging::info("Summary DONE");
}

auto Particle3D::update_direct_force(Eigen::Vector3d force) -> void {
    this->m_direct_force = force;
}

auto Particle3D::calc_direct_force_with_part(const Particle3D &other_part) -> Eigen::Vector3d {
    // (dhub) Source:
    // https://www.physicsforums.com/threads/solving-n-body-simulation-problems-with-gravitational-equations.455058/
    auto dist_part = this->m_position - other_part.m_position;
    auto dist_norm =
        std::pow(dist_part.dot(dist_part) + std::pow(Particle3D::s_softening, 2), 3 / 2);

    return -other_part.km_non_dim_mass * dist_part / dist_norm;
}
