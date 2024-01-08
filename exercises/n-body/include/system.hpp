#ifndef COMPASTRO_SYSTEM_H_
#define COMPASTRO_SYSTEM_H_

// #include "types.hpp"
#include "particle.hpp"
#include "shell.hpp"

class System {
  public:
    PartVec m_particles;
    double m_total_mass = 0.;
    double m_half_mass = 0.;
    double m_scale_length = 0.;
    double m_min_rad = 0.;
    double m_max_rad = 0.;

    System();
    System(System &&) = default;
    System(const System &) = default;
    System &operator=(System &&) = default;
    System &operator=(const System &) = default;
    ~System();

    auto transform_vectors()
        -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>;

    /// @brief Return the particle that is the furthest away
    /// @note In a running system, the distances need to be calculated at each step
    auto get_max_distance() -> Particle3D;

    auto calc_total_mass() const -> double;
    auto update_total_mass() -> void;

    auto calc_half_mass(const ShellVec &shells) const -> double;
    auto update_half_mass(const ShellVec &shells) -> void;

    auto calc_scale_length() const -> double;
    auto update_scale_length() -> void;

    auto update_min_rad(const double rad) -> void;
    auto update_max_rad(const double rad) -> void;

    auto density_hernquist(double rad) -> double;

  private:
};

#endif // ! COMPASTRO_SYSTEM_H_
