#include "system.hpp"
#include "data.hpp"
#include "logging.hpp"
#include "particle.hpp"

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <future>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>

System::System(const std::string_view &path_name) {
    auto particles_opt = Data::read_data(path_name);
    if (not particles_opt.has_value()) {
        Logging::err("Error while reading file: {}", path_name);
        std::exit(-1);
    };
    m_particles = particles_opt.value();
}

// TODO: (aver) Empty constructor needed for later global initialization
System::System() {}
System::~System() {}

auto System::precalc_consts() -> void {
    for (const auto &part : m_particles) {
        update_min_rad(part.m_distance);
        update_max_rad(part.m_distance);
        m_total_mass += Particle3D::km_non_dim_mass;
    }
    // WARN: (aver) We need min rad to be larger than 0, otherwise many following calculations have
    // a divide by 0!
    assert(m_min_rad > 0.);
}

auto System::transform_vectors()
    -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    // TODO: (aver) convert to openmp
    //
    auto future_x = std::async(std::launch::async, [&]() {
        auto transformed_x = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.m_position.x();
                             });
        std::vector<double> transformed_vector_x(transformed_x.begin(), transformed_x.end());
        return transformed_vector_x;
    });
    auto future_y = std::async(std::launch::async, [&]() {
        auto transformed_y = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.m_position.y();
                             });
        std::vector<double> transformed_vector_y(transformed_y.begin(), transformed_y.end());
        return transformed_vector_y;
    });
    auto future_z = std::async(std::launch::async, [&]() {
        auto transformed_z = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.m_position.z();
                             });
        std::vector<double> transformed_vector_z(transformed_z.begin(), transformed_z.end());
        return transformed_vector_z;
    });

    auto x = future_x.get();
    auto y = future_y.get();
    auto z = future_z.get();

    return {x, y, z};
}
// TODO: (aver) add min function
auto System::get_max_distance() -> Particle3D {
    auto result_it =
        std::max_element(m_particles.begin(), m_particles.end(), [](Particle3D a, Particle3D b) {
            return a.m_distance < b.m_distance;
        });
    // result_it->print_summary();
    return *result_it;
}
auto System::calc_total_mass() const -> double {

    auto total_mass = std::accumulate(
        m_particles.begin(), m_particles.end(), 0., [&](double sum, const Particle3D &part) {
            // return sum + part.mass;
            return sum + Particle3D::km_non_dim_mass;
        });
    Logging::info("Total mass of system: {}", m_total_mass);

    return total_mass;
}
auto System::update_total_mass() -> void { m_total_mass = calc_total_mass(); }

auto System::calc_half_mass_radius(const ShellVec &shells) const -> double {

    auto temp_mass = 0.;
    for (const auto &shell : shells) {
        if (temp_mass >= m_total_mass * 0.5) {
            return shell.m_upper;
        }
        temp_mass += shell.m_mass;
    }
    Logging::err("No half mass found... Exiting");
    std::exit(-1);
}

auto System::update_half_mass_radius(const ShellVec &shells) -> void {
    m_half_mass_rad = calc_half_mass_radius(shells);
}

auto System::calc_scale_length() const -> double { return m_half_mass_rad / (1 + std::sqrt(2)); }
auto System::update_scale_length() -> void { m_scale_length = calc_scale_length(); }

auto System::update_min_rad(const double rad) -> void { m_min_rad = std::min(m_min_rad, rad); }
auto System::update_max_rad(const double rad) -> void { m_max_rad = std::max(m_max_rad, rad); }

auto System::get_constrained_shell_mass(const double lower_rad, const double upper_rad) const
    -> double {
    return std::accumulate(
        m_particles.begin(), m_particles.end(), 0., [&](double sum, const Particle3D &part) {
            if (part.m_distance >= lower_rad and part.m_distance <= upper_rad) {
                // return sum + part.mass;
                return sum + Particle3D::km_non_dim_mass;
            }
            return sum + 0.;
        });
}

auto System::density_hernquist(const double rad) const -> double {
    return (m_total_mass * m_scale_length) /
           (2 * std::numbers::pi * rad * std::pow(rad + m_scale_length, 3));
}
auto System::newton_force(const double rad) const -> double {
    return -m_total_mass * Particle3D::km_non_dim_mass / (std::pow(rad + m_scale_length, 2));
}

auto System::calc_direct_force() -> void {
#pragma omp parallel for
    for (uint i = 0; i < m_particles.size(); i++) {
        Eigen::Vector3d sum_force_inter_part({0, 0, 0});
        for (uint j = 0; j < m_particles.size(); j++) {
            if (i != j) {
                sum_force_inter_part += m_particles[i].calc_direct_force_with_part(m_particles[j]);
            }
        }
#pragma omp critical
        m_particles[i].update_direct_force(sum_force_inter_part);
    }
}
auto System::convert_lin_to_log(const int no_bins, const double val) const -> double {
    return m_min_rad * std::pow(m_max_rad / m_min_rad, val / no_bins);
}

auto System::fit_log_to_plot(const double val) -> double {
    return val + std::numeric_limits<double>::epsilon();
}
