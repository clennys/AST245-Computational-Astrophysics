#include "shell.hpp"

Shell::Shell(const uint idx, const double lower_inc, const double upper)
    : m_index(idx), m_lower_inc(lower_inc), m_upper(upper) {}

auto Shell::shell_int_size() const -> int { return static_cast<int>(m_particles.size()); }
auto Shell::calc_volume() const -> double {
    return k_shell_vol_pref * (std::pow(m_upper, 3) - std::pow(m_lower_inc, 3));
}

auto Shell::update_volume() -> void { m_volume = calc_volume(); }

auto Shell::calc_density() const -> double { return m_mass / m_volume; }

auto Shell::update_density() -> void { m_density = calc_density(); }
