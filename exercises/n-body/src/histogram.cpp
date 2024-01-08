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

Histogram::Histogram(const uint no_bins, const double radius, System &p_system) {

    // WARN: (aver) Addded one to radious to securely place all particles into shells later on
    auto bin_radius = (radius + 1) / static_cast<int>(no_bins);
    auto lower_rad = 0.;
    auto upper_rad = bin_radius;

    // setup bins/shells
    for (uint i = 0; i < no_bins; i++) {
        // Logging::dbg(lower_rad);
        // Logging::dbg(upper_rad);

        // m_shells.push_back(Shell(lower_rad, upper_rad));
        m_shells.emplace_back(Shell(i, lower_rad, upper_rad));
        lower_rad += bin_radius;
        upper_rad += bin_radius;
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
        p_system.m_total_mass += part.mass;
        p_system.update_min_rad(part.distance);
        // this is calculated before calling the histogram constructor
        p_system.update_max_rad(part.distance);

        auto shell_it = std::find_if(m_shells.begin(), m_shells.end(), [&part](const Shell &shell) {
            // Logging::dbg(std::format("Working on shell lower: {}", shell.m_lower_inc));
            // Logging::dbg(std::format("Working on shell upper: {}", shell.m_upper));
            // Logging::dbg(std::format("Working on particle with distance: {}", part.distance));
            return part.distance < shell.m_upper && shell.m_lower_inc >= part.distance;
        });

        if (shell_it != m_shells.end()) {
            // Logging::dbg(std::format("Working on shell: {}", it->m_index));
            shell_it->m_particles.emplace_back(part);
            shell_it->m_mass += part.mass;
        } else {
            // handle the case, where a particle is not placed into any of the shells...
            Logging::err(
                "Particle with distance: {}, was not placed into a shell.\n\tShell bound: [{}, {})",
                part.distance,
                shell_it->m_lower_inc,
                shell_it->m_upper);

            // std::exit(-1);
        }
    }
#endif

    Logging::info("...Particles emplaced!");
}

Histogram::~Histogram() {}
