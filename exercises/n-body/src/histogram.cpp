#include "histogram.hpp"
#include "logging.hpp"
#include "shell.hpp"
#include "types.hpp"

#include <algorithm>
#include <execution>
#include <format>
#include <mutex>

Histogram::Histogram(const uint no_bins, const double radius, const PartVec &particles) {
    auto bin_radius = radius / static_cast<int>(no_bins);
    auto lower_rad = 0.;
    auto upper_rad = bin_radius;

    std::mutex shell_mutex;

    // setup bins/shells
    for (uint i = 0; i < no_bins; i++) {
        // Logging::dbg(lower_rad);
        // Logging::dbg(upper_rad);

        // WARN: (aver) don't do this, ends in segfault, because of double rounding operation
        // if (i == no_bins - 1) {
        //     // in case we are at the last bin, we round up the upper bound
        //     upper_rad = std::ceil(upper_rad);
        // }

        // m_shells.push_back(Shell(lower_rad, upper_rad));
        m_shells.emplace_back(Shell(i, lower_rad, upper_rad));
        lower_rad += bin_radius;
        upper_rad += bin_radius;
    }

    Logging::info("Sorting Particles into shells...");

#if 0
    std::for_each(
        std::execution::par,
        particles.begin(),
        particles.end(),
        [this, &shell_mutex](const Particle3D &part) {
            auto it = std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
                return part.distance < shell.m_upper && shell.m_lower_inc >= part.distance;
            });

            std::lock_guard<std::mutex> guard(shell_mutex); // Protect m_shells
            if (it != m_shells.end()) {
                it->m_particles.emplace_back(part);
            } else {
                // handle the case, where a particle is not placed into any of the shells...
                Logging::err(std::format("Particle with distance: {}, was not placed into a "
                                         "shell.\n\tShell bound: [{}, {})",
                                         part.distance,
                                         it->m_lower_inc,
                                         it->m_upper));
                // std::exit(-1);
            }
        });
#else
    for (const auto &part : particles) {
        auto it = std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
            // Logging::dbg(std::format("Working on shell lower: {}", shell.m_lower_inc));
            // Logging::dbg(std::format("Working on shell upper: {}", shell.m_upper));
            // Logging::dbg(std::format("Working on particle with distance: {}", part.distance));
            return part.distance < shell.m_upper && shell.m_lower_inc >= part.distance;
        });

        if (it != m_shells.end()) {
            // Logging::dbg(std::format("Working on shell: {}", it->m_index));
            it->m_particles.emplace_back(part);
        } else {
            // handle the case, where a particle is not placed into any of the shells...
            Logging::err(std::format(
                "Particle with distance: {}, was not placed into a shell.\n\tShell bound: [{}, {})",
                part.distance,
                it->m_lower_inc,
                it->m_upper));

            // std::exit(-1);
        }
    }
#endif

    Logging::info("DONE");
}

Histogram::~Histogram() {}
