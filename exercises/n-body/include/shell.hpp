#ifndef SHELL_H_
#define SHELL_H_
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

    /// Hold particles that are in the shell
    PartVec m_particles;

    explicit Shell(const uint idx, const double lower_inc, const double upper);
    Shell(Shell &&) = default;
    Shell(const Shell &) = default;
    Shell &operator=(Shell &&) = default;
    Shell &operator=(const Shell &) = default;
    ~Shell();

  private:
};

#endif // ! SHELL_H_
