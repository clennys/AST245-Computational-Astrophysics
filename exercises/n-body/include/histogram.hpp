#ifndef HISTOGRAM_H_
#define HISTOGRAM_H_

#include "shell.hpp"
#include "system.hpp"

/// Creates a histogram with `no_bins` shells. If `do_log` is `true`, then the size of the bins
/// is logarithmically placed.
class Histogram {
  public:
    ShellVec m_shells;

    explicit Histogram(const int no_bins, const System &p_system, bool do_log = false);
    Histogram(Histogram &&) = default;
    Histogram(const Histogram &) = default;
    Histogram &operator=(Histogram &&) = default;
    Histogram &operator=(const Histogram &) = default;
    ~Histogram() = default;

  private:
};

#endif // ! HISTOGRAM_H_
