#ifndef COMPASTRO_SHELL_H_
#define COMPASTRO_SHELL_H_

#include "particle.hpp"
#include "types.hpp"

class Shell {
  public:
    /// index of current shell
    uint m_index;
    /// Lower range, inclusive
    double m_lower_inc;
    /// Upper range, non-inclusive
    double m_upper;
    /// Hold the mass of all the particles in the shell
    double m_mass = 0.;
    double m_volume;
    double m_density;

    /// Hold particles that are in the shell
    PartVec m_particles;

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

  private:
    // precalculated (by the compiler!) constant for shell volume
    static constexpr auto k_shell_vol_pref = 4. / 3. * std::numbers::pi;
};

using ShellVec = std::vector<Shell>;

#endif // ! COMPASTRO_SHELL_H_
