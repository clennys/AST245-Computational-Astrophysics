#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "shell.hpp"
#include "system.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstdlib>
#include <execution>
#include <format>
#include <mutex>

Histogram::Histogram(const int no_bins, const double radius, System &p_system, bool do_log) {

    // WARN: (aver) Addded one to radious to securely place all particles into shells later on
    auto bin_radius = (radius + 1) / static_cast<int>(no_bins);
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

        m_shells.emplace_back(Shell(i, lower_rad, upper_rad));
    }

    Logging::info("Sorting Particles into shells...");

#if 0
    std::mutex shell_mutex;

    std::for_each(
        std::execution::par,
        p_system.m_particles.begin(),
        p_system.m_particles.end(),
        [this, &shell_mutex, &p_system](const Particle3D &part) {
            // use this loop to get total mass
            p_system.m_total_mass += part.mass;
            p_system.update_min_rad(part.distance);

            auto shell_it =
                std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
                    return part.distance < shell.m_upper && shell.m_lower_inc >= part.distance;
                });

            std::lock_guard<std::mutex> guard(shell_mutex); // Protect m_shells
            // use this loop to get total mass
            p_system.m_total_mass += part.mass;

            if (shell_it != m_shells.end()) {
                shell_it->m_particles.emplace_back(part);
                shell_it->m_mass += part.mass;
            } else {
                // handle the case, where a particle is not placed into any of the shells...
                Logging::err(std::format("Particle with distance: {}, was not placed into a "
                                         "shell.\n\tShell bound: [{}, {})",
                                         part.distance,
                                         shell_it->m_lower_inc,
                                         shell_it->m_upper));
            }
        });
#else
    for (const auto &part : p_system.m_particles) {

        // use this loop to get total mass
        // p_system.m_total_mass += part.mass;
        p_system.m_total_mass += Particle3D::km_non_dim_mass;

        p_system.update_min_rad(part.m_distance);
        // this is calculated before calling the histogram constructor
        p_system.update_max_rad(part.m_distance);

        auto shell_it = std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
            // Logging::dbg(std::format("Working on shell lower: {}", shell.m_lower_inc));
            // Logging::dbg(std::format("Working on shell upper: {}", shell.m_upper));
            // Logging::dbg(std::format("Working on particle with distance: {}", part.distance));
            return part.m_distance < shell.m_upper && shell.m_lower_inc >= part.m_distance;
        });

        if (shell_it != m_shells.end()) {
            // Logging::dbg(std::format("Working on shell: {}", it->m_index));
            shell_it->m_particles.emplace_back(part);
            // shell_it->m_mass += part.mass;
            shell_it->m_mass += Particle3D::km_non_dim_mass;
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
#endif

    Logging::info("...Particles emplaced!");
}

Histogram::~Histogram() {}
