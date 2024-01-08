#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <Eigen/Dense>

class Particle3D {
  public:
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double softening;
    double potential;
    double distance;
    double mass;

    /// Calculate the Norm of the position vector and set the value
    auto calc_orign_distance() -> void;

    /// Print a summary of the particle properties
    auto print_summary() const -> void;
};

using PartVec = std::vector<Particle3D>;
#endif // ! PARTICLES_H_
