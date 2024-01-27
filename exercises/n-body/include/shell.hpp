#ifndef COMPASTRO_SHELL_H_
#define COMPASTRO_SHELL_H_

#include "particle.hpp"
#include "types.hpp"

class Shell {
  public:
    /// index of current shell
    uint m_index = 0;
    /// Lower range, inclusive
    double m_lower_inc = 0;
    /// Upper range, non-inclusive
    double m_upper = 0;
    /// Hold the mass of all the particles in the shell
    double m_mass = 0.;
    double m_volume = 0;
    double m_density = 0;

    /// Hold particles that are in the shell
    PartVec m_particles = {};

    explicit Shell(const uint idx, const double lower_inc, const double upper);
    Shell(Shell &&) = default;
    Shell(const Shell &) = default;
    Shell &operator=(Shell &&) = default;
    Shell &operator=(const Shell &) = default;
    ~Shell() = default;

    /// Return the count of particles in the shell
    auto shell_int_size() const -> int;

    auto calc_volume() const -> double;

    auto update_volume() -> void;

    auto calc_density() const -> double;

    auto update_density() -> void;

    /// Return the average force calculated bu direct summation in the shell, projected on the
    /// center of the spherical distribution
    auto get_avg_direct_force() const -> double;
    auto get_avg_tree_force() const -> double;

  private:
    // precalculated (by the compiler!) constant for shell volume
    static constexpr auto k_shell_vol_pref = 4. / 3. * std::numbers::pi;
};

using ShellVec = std::vector<Shell>;

#endif // ! COMPASTRO_SHELL_H_
