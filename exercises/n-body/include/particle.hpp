#ifndef PARTICLES_H_
#define PARTICLES_H_

#include <Eigen/Dense>

class Particle3D {
  public:
    double mass;
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double softening;
    double potential;
};

#endif // ! PARTICLES_H_
