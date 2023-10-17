#pragma once
#include <Eigen/Dense>

class Particle {
  public:
    Eigen::Vector2d position;
    Eigen::Vector2d velocity;
    const double k_mass;
    // TODO: (aver) Consider them making private
    double energy;
    double eccentricity;
    double angular_momentum;
    double time;

    Particle(Eigen::Vector2d pos, Eigen::Vector2d vel, double mass);
    Particle(Particle &&) = default;
    Particle(const Particle &) = default;
    Particle &operator=(Particle &&) = default;
    Particle &operator=(const Particle &);
    ~Particle();

    auto compute_energy() -> void;
    auto compute_angular_momentum() -> void;
    auto compute_eccentricity() -> void;

    // TODO: (aver) implement operator for printing.
    // auto operator<<(std::ostream &os);

  private:
};
