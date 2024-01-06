#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <Eigen/Dense>
#include <vector>


class Particle3D {
  public:
    double mass;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double softening;
    double potential;
    double distance;

    /// @brief calculate the Norm of the position vector
    auto calc_orign_distance() -> void;
		auto print_summary() -> void;
};

/// Utility functions to do work with a `std::vector<Particle3D>`
namespace Particles {

/// @brief Return the particle that is the furthest away
/// @note In a running system, the distances need to be calculated at each step
auto get_max_distance(const std::vector<Particle3D> &particles) -> Particle3D;
} // namespace Particles


#endif // ! PARTICLES_H_
