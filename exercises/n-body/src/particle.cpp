#include "particle.hpp"
#include "logging.hpp"
#include "system.hpp"

auto Particle3D::calc_origin_dist() -> double { return this->m_position.norm(); }

auto Particle3D::update_origin_dist() -> void { this->m_distance = calc_origin_dist(); }

auto Particle3D::print_summary() const -> void {
    Logging::info("Printing particle summary:");
    std::cout << "Position =\n" << this->m_position << std::endl;
    std::cout << "Velocity =\n" << this->m_velocity << std::endl;
    std::cout << "Potential =\n" << this->m_potential << std::endl;
    std::cout << "Distance =\n" << this->m_distance << std::endl;
    std::cout << "Force =\n" << this->m_direct_force << std::endl;
}

auto Particle3D::update_direct_force(const Eigen::Vector3d force) -> void {
    this->m_direct_force = force;
}

auto Particle3D::calc_direct_force_with_part(const Particle3D &other_part) const
    -> Eigen::Vector3d {
    // (dhub) Source:
    // https://www.physicsforums.com/threads/solving-n-body-simulation-problems-with-gravitational-equations.455058/

    const auto diff_part = this->m_position - other_part.m_position;

    // add softening to (^r_ij)^2
    const auto dist_norm =
        std::sqrt(diff_part.squaredNorm() + (System::s_softening * System::s_softening));

    const auto force_magn =
        (this->m_mass * other_part.m_mass) / (dist_norm * dist_norm * dist_norm);

    // assumption: G=1, and Mass=1, otherwise the mass of Particle_i would need to be multiplied
    // here as well
    return -force_magn * diff_part;
}

auto Particle3D::historize_part_state() -> void { m_pos_history.push_back(m_position); }
