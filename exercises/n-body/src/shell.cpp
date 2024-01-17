#include "shell.hpp"

Shell::Shell(const uint idx, const double lower_inc, const double upper)
    : m_index(idx), m_lower_inc(lower_inc), m_upper(upper) {}

Shell::~Shell() {}

auto Shell::shell_int_size() const -> int { return static_cast<int>(m_particles.size()); }
