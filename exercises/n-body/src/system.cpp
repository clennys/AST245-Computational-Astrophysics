#include "system.hpp"
#include "logging.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <future>
#include <numbers>
#include <numeric>
#include <ranges>
#include <tuple>

System::System() {}

System::~System() {}

auto System::transform_vectors()
    -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    // TODO: (aver) convert to openmp
    //
    auto future_x = std::async(std::launch::async, [&]() {
        auto transformed_x = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.x();
                             });
        std::vector<double> transformed_vector_x(transformed_x.begin(), transformed_x.end());
        return transformed_vector_x;
    });
    auto future_y = std::async(std::launch::async, [&]() {
        auto transformed_y = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.y();
                             });
        std::vector<double> transformed_vector_y(transformed_y.begin(), transformed_y.end());
        return transformed_vector_y;
    });
    auto future_z = std::async(std::launch::async, [&]() {
        auto transformed_z = m_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.z();
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
    auto result =
        std::max_element(m_particles.begin(), m_particles.end(), [](Particle3D a, Particle3D b) {
            return a.distance < b.distance;
        });

    result->print_summary();

    return *result;
}
auto System::calc_total_mass() const -> double {

    auto total_mass = std::accumulate(
        m_particles.begin(), m_particles.end(), 0., [&](double sum, const Particle3D &part) {
            return sum + part.mass;
        });

    //     auto total_mass = 0.;
    // #pragma omp parallel for
    //     for (const auto &part : m_particles) {
    // #pragma omp atomic
    //         total_mass += part.mass;
    //     }
    Logging::info(std::format("Total mass of system: {}", m_total_mass));

    return total_mass;
}
auto System::update_total_mass() -> void { m_total_mass = calc_total_mass(); }

auto System::calc_half_mass(const ShellVec &shells) const -> double {

    auto temp_mass = 0.;
    for (const auto &shell : shells) {
        if (temp_mass >= m_total_mass * 0.5) {
            return temp_mass;
        }
        temp_mass += shell.m_mass;
    }
    Logging::err("No half mass found... Exiting");
    std::exit(-1);
}

auto System::update_half_mass(const ShellVec &shells) -> void {
    m_half_mass = calc_half_mass(shells);
}

auto System::calc_scale_length() const -> double { return m_half_mass / (1 / std::sqrt(2)); }
auto System::update_scale_length() -> void { m_scale_length = calc_scale_length(); }

auto System::update_min_rad(const double rad) -> void { m_min_rad = std::min(m_min_rad, rad); }
auto System::update_max_rad(const double rad) -> void { m_max_rad = std::max(m_max_rad, rad); }

auto System::density_hernquist(double rad) -> double {
    return (m_total_mass / (2 * std::numbers::pi)) * (m_scale_length / rad) *
           (1 / std::pow(rad + m_scale_length, 3));
}
