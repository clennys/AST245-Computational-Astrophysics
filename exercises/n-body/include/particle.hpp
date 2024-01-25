#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <Eigen/Dense>

class Particle3D {
  public:
    //=============================================================================================
    // Regular member variables
    //=============================================================================================
    uint m_id;
    Eigen::Vector3d m_position;
    Eigen::Vector3d m_velocity;
    Eigen::Vector3d m_direct_force;
    Eigen::Vector3d m_tree_force;
    double m_potential = 0;
    double m_distance = 0;
    double m_mass = 0;
    // TODO: (dhub) Change input, consider removing
    // double m_softening=-1;
    double m_treecode_potential = 0;

    /// Calculate the Norm of the position vector and return it
    auto calc_origin_dist() -> double;
    /// Internally calculate and set member variable
    auto update_origin_dist() -> void;

    /// Print a summary of the particle properties
    auto print_summary() const -> void;

    /// Calculate direct force with another particle
    auto calc_direct_force_with_part(const Particle3D &other_part) const -> Eigen::Vector3d;

    /// Update direct force on particle
    auto update_direct_force(Eigen::Vector3d force) -> void;

    auto norm_force() -> double;
};

using PartVec = std::vector<Particle3D>;
#endif // ! PARTICLES_H_
