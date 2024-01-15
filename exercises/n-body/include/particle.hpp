#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <Eigen/Dense>

class Particle3D {
  public:
    static constexpr double km_non_dim_mass = 1.;
    // static constexpr double km_non_dim_mass = 92.4259;
    static double s_softening;

    Eigen::Vector3d m_position;
    Eigen::Vector3d m_velocity;
    Eigen::Vector3d m_direct_force;
    double m_potential;
    double m_distance;
    double m_mass;
    // TODO: (dhub) Change input, consider removing
    // double m_softening;

    /// Calculate the Norm of the position vector and return it
    auto calc_origin_dist() -> double;
    /// Internally calculate and set member variable
    auto update_origin_dist() -> void;

    /// Print a summary of the particle properties
    auto print_summary() const -> void;

    /// Calculate direct force with another particle
    auto calc_direct_force_with_part(const Particle3D &other_part) -> Eigen::Vector3d;

    /// Update direct force on particle
    auto update_direct_force(Eigen::Vector3d force) -> void;

    auto norm_force() -> double;
};

using PartVec = std::vector<Particle3D>;
#endif // ! PARTICLES_H_
