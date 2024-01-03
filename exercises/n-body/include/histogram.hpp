#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "shell.hpp"
#include "types.hpp"
#include <vector>

class Histogram {
  public:
    explicit Histogram(const int no_bins, const PartVec &particles);
    Histogram(Histogram &&) = default;
    Histogram(const Histogram &) = default;
    Histogram &operator=(Histogram &&) = default;
    Histogram &operator=(const Histogram &) = default;
    ~Histogram();

  private:
    std::vector<Shell> m_shells;
};

#endif // ! HISTOGRAM_H_
