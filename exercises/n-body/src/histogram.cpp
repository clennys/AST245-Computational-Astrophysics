#include "histogram.hpp"
#include "shell.hpp"

Histogram::Histogram(const int no_bins, const PartVec &particles) {
    for (size_t i = 0; i < no_bins; i++) {
        m_shells.push_back(Shell());
    }
}

Histogram::~Histogram() {}
