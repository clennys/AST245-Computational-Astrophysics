#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "shell.hpp"
#include "system.hpp"
#include "types.hpp"

class Histogram {
  public:
    ShellVec m_shells;

    explicit Histogram(const uint no_bins, const double radius, System &p_system);
    Histogram(Histogram &&) = default;
    Histogram(const Histogram &) = default;
    Histogram &operator=(Histogram &&) = default;
    Histogram &operator=(const Histogram &) = default;
    ~Histogram();

  private:
};

#endif // ! HISTOGRAM_H_
