#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "shell.hpp"
#include "system.hpp"
#include "types.hpp"

#include <algorithm>
#include <format>
#include <omp.h>

Histogram::Histogram(const int no_bins, const System &p_system, bool do_log) {
    // Logging::info("Creating Histogram with {} shells...", no_bins);

    // WARN: (aver) Addded one to radious to securely place all particles into shells later on
    const auto bin_radius = (p_system.m_max_rad + 1) / no_bins;

    double lower_rad;
    double upper_rad;

    // setup bins/shells
    for (int i = 0; i <= no_bins; i++) {
        if (do_log) {
            lower_rad = p_system.convert_lin_to_log(no_bins, i);
            upper_rad = p_system.convert_lin_to_log(no_bins, i + 1);
        } else {
            lower_rad = i * bin_radius;
            upper_rad = (i + 1) * bin_radius;
        }

        m_shells.emplace_back(Shell(static_cast<uint>(i), lower_rad, upper_rad));
    }

    for (const auto &part : p_system.m_particles) {
        auto shell_it = std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
            return part.m_distance >= shell.m_lower_inc and part.m_distance < shell.m_upper;
        });

        if (shell_it != m_shells.end()) {
            // Logging::dbg(std::format("Working on shell: {}", it->m_index));
            shell_it->m_particles.emplace_back(part);
            shell_it->m_mass += part.m_mass;
        } else {
            // handle the case, where a particle is not placed into any of the shells...
            Logging::err(
                "Particle with distance: {}, was not placed into a shell.\n\tShell bound: [{}, {})",
                part.m_distance,
                shell_it->m_lower_inc,
                shell_it->m_upper);

            // std::exit(-1);
        }
    }

    // We loop now through all shells und calculate their volumes and densities
#pragma omp parallel for
    for (auto &shell : m_shells) {
        shell.update_volume();
        shell.update_density();
    }
}
