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

// set default softening
double System::s_softening = 0.;

System::System(const std::string_view &path_name) {
    auto particles_opt = Data::read_data(path_name);
    if (not particles_opt.has_value()) {
        Logging::err("Error while reading file: {}", path_name);
        std::exit(-1);
    };
    m_particles = particles_opt.value();
}

auto System::system_int_size() const -> int { return static_cast<int>(m_particles.size()); }

auto System::convert_lin_to_log(const int no_bins, const double val) const -> double {
    return m_min_rad * std::pow(m_max_rad / m_min_rad, val / no_bins);
}

auto System::fit_log_to_plot(const double val) -> double {
    return val + std::numeric_limits<double>::epsilon();
}

auto System::precalc_consts() -> void {
    for (const auto &part : m_particles) {
        update_min_rad(part.m_distance);
        update_max_rad(part.m_distance);
        m_total_mass += System::k_non_dim_mass;
    }
    // WARN: (aver) We need min rad to be larger than 0, otherwise many following calculations have
    // a divide by 0!
    assert(m_min_rad > 0.);
}

auto System::precalc_mean_inter_part_dist() -> double {
    auto mean_dist = 0.;

#pragma omp parallel for reduction(+ : mean_dist)
    for (uint64_t i = 0; i < m_particles.size(); ++i) {
#pragma omp parallel for reduction(+ : mean_dist)
        for (uint64_t j = 0; j < m_particles.size(); ++j) {
            if (i == j)
                continue;
            mean_dist += (m_particles[i].m_position - m_particles[j].m_position).norm();
        }
    }
    // WARN: (aver) casting to double results in an illegal instruction, so bear with this
    // 64bit needed to hold the value, which is larger than 2^32
    // mean_dist /= static_cast<int64_t>(m_particles.size() * (m_particles.size() - 1));
    mean_dist *= 2. / static_cast<double>(
                          static_cast<int64_t>(m_particles.size() * (m_particles.size() - 1)));

    return mean_dist;
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

            (void)part;
            return sum + System::k_non_dim_mass;
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

auto System::density_hernquist(const double rad) const -> double {
    return (m_total_mass * m_scale_length) /
           (2 * std::numbers::pi * rad * std::pow(rad + m_scale_length, 3));
}
auto System::newton_force(const double rad) const -> double {
    // applying Newtons second theoreom on spherical systems and the M(r) function of the Hernquist
    // paper, we get this calculation, while rad^2 can be reduced

    // auto M_r = m_total_mass * (rad * rad) / ((rad + m_scale_length) * (rad + m_scale_length));
    // return -M_r / (rad * rad);

    return -m_total_mass * System::k_non_dim_mass /
           ((rad + m_scale_length) * (rad + m_scale_length));
}

auto System::calc_direct_force() -> void {
#if 1
#pragma omp parallel for
    for (uint64_t i = 0; i < m_particles.size(); ++i) {
        // reset on each `i` change
        auto sum_force_inter_part = Eigen::Vector3d({0, 0, 0});

        for (uint64_t j = 0; j < m_particles.size(); ++j) {
            if (i == j)
                continue;

            auto val = m_particles[i].calc_direct_force_with_part(m_particles[j]);
            sum_force_inter_part += val;
        }
#pragma omp critical
        m_particles[i].update_direct_force(sum_force_inter_part);
    }
#else
    // No perf gain in parallelizing this ...
    for (uint i = 0; i < m_particles.size(); ++i) {
        for (uint j = i + 1; j < m_particles.size(); ++j) {
            const auto force_vec = m_particles[i].calc_direct_force_with_part(m_particles[j]);
            m_particles[i].m_direct_force += force_vec;
            m_particles[j].m_direct_force -= force_vec;
        }
    }
#endif
}

auto System::solver_do_step(const double delta_time) -> void {
    // Logging::info("Stepping forward with dt: {}", delta_time);

#pragma omp parallel for
    for (auto &part : m_particles) {
        const auto velocity_mid =
            part.m_velocity + (part.m_direct_force / System::k_non_dim_mass) * delta_time / 2;

        part.m_position += velocity_mid * delta_time;
        auto force_new = Eigen::Vector3d({0, 0, 0});
        // update forces at new position

#if 1
#pragma omp parallel private(part)
        {
            auto local_force = Eigen::Vector3d({0, 0, 0});

#pragma omp for
            for (size_t i = 0; i < m_particles.size(); ++i) {
                auto other_part = m_particles[i];
                if (other_part.m_id == part.m_id)
                    continue;

                local_force += part.calc_direct_force_with_part(other_part);
            }
#pragma omp critical
            force_new += local_force;
        }
#else
#pragma omp parallel for private(part)
        for (size_t i = 0; i < m_particles.size(); i++) {
            auto other_part = m_particles[i];
            if (other_part.m_id == part.m_id)
                continue;
#pragma omp critical
            force_new += part.calc_direct_force_with_part(other_part);
        }
#endif

        // update the force
        part.update_direct_force(force_new);

        // set new velocity, completing leap-frog
        part.m_velocity =
            velocity_mid + (part.m_direct_force / System::k_non_dim_mass) * delta_time / 2;
    }
}

auto System::calc_relaxation() const -> double {
    // NOTE: (dhub) Assume G=1
    double nr_part = this->system_int_size();
    double circular_velocity = std::sqrt(m_total_mass * 0.5 / m_half_mass_rad);
    double time_cross = m_half_mass_rad / circular_velocity;
    return nr_part / (8 * std::log(nr_part)) * time_cross;
}

auto System::update_relaxation() -> void { m_relaxation = calc_relaxation(); }
