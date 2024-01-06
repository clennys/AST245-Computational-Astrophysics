#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "shell.hpp"
#include "types.hpp"
#include <vector>

class Histogram {
  public:
    std::vector<Shell> m_shells;

    explicit Histogram(const uint no_bins, const double radius, const PartVec &particles);
    Histogram(Histogram &&) = default;
    Histogram(const Histogram &) = default;
    Histogram &operator=(Histogram &&) = default;
    Histogram &operator=(const Histogram &) = default;
    ~Histogram();

  private:
};

#endif // ! HISTOGRAM_H_
