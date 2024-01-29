#include "system.hpp"
#include "data.hpp"
#include "histogram.hpp"
#include "logging.hpp"
#include "mgl2/mgl.h"
#include "node.hpp"
#include "particle.hpp"
#include <sstream>

#include "Eigen/Eigen"

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
double System::s_softening_length = 0.;

auto System::init_system(const std::string_view &path_name) -> void {
    Logging::info("");
    Logging::info("==============================================================================");
    Logging::info("Initializing System with file: {}", path_name);

    auto particles_opt = Data::read_data(path_name);
    if (not particles_opt.has_value()) {
        Logging::err("Error while reading file: {}", path_name);
        std::exit(-1);
    };
    this->ic_particles = particles_opt.value();
    this->m_particles = particles_opt.value();

    this->precalc_consts();

    // use this shell to calculcate half_mass_radius, and scale_length
    Histogram hist(100'000, *this);

    this->update_half_mass_radius(hist.m_shells);
    this->update_scale_length();
    this->update_times();

    // Relaxation under constraints is not needed
    // this->update_relaxation();

    /* In collisionless simulations the choice of softening kernel and softening
    length is a compromise between maximizing the smoothness of the force-
    field and minimizing the degradation of the spatial resolution caused by the
    softening. Despite its popularity, the choice (2.226) of S is not optimal in
    this sense, because the density of a Plummer sphere falls off too slowly with
    radius */

    // System::s_softening_length =
    //     -1 / (std::sqrt(this->m_max_rad * this->m_max_rad +
    //                     System::k_mean_inter_dist * System::k_mean_inter_dist));

    // System::s_softening_length =
    //     -(this->m_max_rad * this->m_max_rad +
    //       (3. / 2.) * System::k_mean_inter_dist * System::k_mean_inter_dist) /
    //     std::pow(std::sqrt(this->m_max_rad * this->m_max_rad +
    //                        System::k_mean_inter_dist * System::k_mean_inter_dist),
    //              3);

    System::s_softening_length =
        -(this->m_half_mass_rad * this->m_half_mass_rad +
          (3. / 2.) * System::k_mean_inter_dist * System::k_mean_inter_dist) /
        std::pow(std::sqrt(this->m_half_mass_rad * this->m_half_mass_rad +
                           System::k_mean_inter_dist * System::k_mean_inter_dist),
                 3);

    System::s_softening = s_softening_length;

    Logging::info("Total mass of system:       {:<12}", m_total_mass);
    Logging::info("Half mass radius of system: {:>12.10f}", m_half_mass_rad);
    Logging::info("Scaling length of system:   {:>12.10f}", m_scale_length);
    Logging::info("Avg Inter-particle dist:    {:>12.10f}", System::k_mean_inter_dist);
    Logging::info("Softening of system:        {:>12.10f}", System::s_softening);
    Logging::info("Crossing time:              {:>12.10f}", m_t_cross);
    Logging::info("Relaxation time:            {:>12.10f}", m_t_relaxation);
    Logging::info("==============================================================================");
    Logging::info("");
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
        m_total_mass += part.m_mass;
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
            return sum + part.m_mass;
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

    return -m_total_mass * System::k_dim_mass / ((rad + m_scale_length) * (rad + m_scale_length));
}

auto System::precalc_direct_initial_force() -> void {
    // TODO: (aver) we could merge this with the solver step, where the initial calculation ignores
    // updating position and velocity
    Logging::info("Calculating initial direct forces...");
    Logging::info("Softening: {}", System::s_softening);
#if 1
    // with omp the 'slower' comparison based loop is much faster with omp (~2s)
#pragma omp parallel for
    for (uint64_t i = 0; i < m_particles.size(); ++i) {
        // reset on each `i` change
        auto sum_force_inter_part = Eigen::Vector3d().setZero();

        for (uint64_t j = 0; j < m_particles.size(); ++j) {
            if (i == j)
                continue;

            auto val = m_particles[i].calc_direct_force_with_part(m_particles[j]);
            sum_force_inter_part += val;
        }
#pragma omp critical
        {
            m_particles[i].update_direct_force(sum_force_inter_part);
            m_particles[i].historize_part_state();
        }
    }
#else
    // No perf gain in parallelizing this ... with omp
    // On a single thread though, proves much more performant (~8 seconds)
    for (uint i = 0; i < m_particles.size(); ++i) {
        for (uint j = i + 1; j < m_particles.size(); ++j) {
            const auto force_vec = m_particles[i].calc_direct_force_with_part(m_particles[j]);
            m_particles[i].m_direct_force += force_vec;
            m_particles[j].m_direct_force -= force_vec;
        }
    }
#endif
}

/// Kick-Drift-Kick
#define KDK

auto System::solver_do_step(const double delta_time) -> void {
#pragma omp parallel for schedule(dynamic)
    for (auto &part : m_particles) {

#ifdef KDK
        // First leap, get mid velocity
        const auto velocity_mid =
            part.m_velocity + (part.m_direct_force / part.m_mass) * delta_time * .5;
        // Update position for force calculation
        part.m_position += velocity_mid * delta_time;
#else
        // r half
        part.m_position = part.m_position + .5 * delta_time * part.m_velocity;
#endif

        auto force_new = Eigen::Vector3d().setZero();
        // update forces at new position
        for (const auto &other : m_particles) {
            if (other.m_id == part.m_id)
                continue;
            force_new += part.calc_direct_force_with_part(other);
        }
        // update the force
        part.update_direct_force(force_new);

#ifdef KDK
        // set new velocity, completing leap-frog
        part.m_velocity = velocity_mid + (part.m_direct_force / part.m_mass) * delta_time * .5;

#else
        // v n+1
        part.m_velocity = part.m_velocity + delta_time * (part.m_direct_force / part.m_mass);
        part.m_position = part.m_position + .5 * delta_time * part.m_velocity;
#endif
        part.historize_part_state();
    }
}

auto System::calc_crossing_time() const -> double {
    double circular_velocity = std::sqrt(m_total_mass * 0.5 / m_half_mass_rad);
    return m_half_mass_rad / circular_velocity;
}

auto System::calc_relaxation_time() const -> double {
    // NOTE: (dhub) Assume G=1
    double nr_part = this->system_int_size();
    // NOTE: (aver) the Mass(R_hm) == half the mass obviously, including the assumption of G==1,
    // m==1
    // Since the half mass radius is fixed, we do not expect change of the relaxation time if the
    // softening were to be increased or decreased. Although numerically, I don't know what the
    // implication are if softening is increased above mean inter particle separation

    // double circular_velocity = std::sqrt(m_total_mass * 0.5 / m_half_mass_rad);
    // double time_cross = m_half_mass_rad / circular_velocity;
    return nr_part / (8 * std::log(nr_part)) * m_t_cross;
}

auto System::update_times() -> void {
    m_t_cross = calc_crossing_time();
    m_t_relaxation = calc_relaxation_time();
}

auto System::calc_overall_bounding_cube() -> BoundingCube {
    // WARN: (aver) Instantiating an Eigen::Tensor mallocs a few times, and apparently it does so
    // asynchronously, leading to an illegal instruction error at cube access
    auto cube = BoundingCube(2, 2, 2);
    auto vec = Eigen::Vector3d().setZero();

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                vec.x() = i == 1 ? m_max_rad : -m_max_rad;
                vec.y() = j == 1 ? m_max_rad : -m_max_rad;
                vec.z() = k == 1 ? m_max_rad : -m_max_rad;
                cube(i, j, k) = vec;
            }
        }
    }

    return cube;
}
auto System::reset_system() -> void { this->m_particles = this->ic_particles; }

auto System::calc_real_relaxation() const -> void {
    Logging::info("==============================================================================");
    Logging::info("Calculating relaxation timescale");

    constexpr auto k_parsec_to_km_factor = 3.0857e+13;
    constexpr auto k_seconds_to_year = 60. * 60. * 24. * 365.25;

    auto rel_hist = Histogram(100'000, *this);
    for (auto &shell : rel_hist.m_shells) {
        shell.m_mass = shell.shell_int_size() * System::k_dim_mass;
        shell.update_density();
    }

    const auto r_hm = this->calc_half_mass_radius(rel_hist.m_shells);
    const auto M = this->system_int_size() * System::k_dim_mass;
    Logging::info("Half Mass Radius at     {} pc", r_hm);

    // WARN: (aver) What to do? Use G=1 or use the real unit of G?
    const auto v_c = std::sqrt(System::k_G * M * 0.5 / r_hm);
    Logging::info("Circular Velocity at    {} km / s", v_c);

    // const auto t_cross = r_hm / v_c;
    // Logging::info("Crossing Timescale at   {} pc / (km/s)", t_cross);
    const auto t_cross_km = (r_hm * k_parsec_to_km_factor) / v_c;
    const auto t_cross_yr = t_cross_km / k_seconds_to_year;
    Logging::info("Crossing Timescale at   {} yr", t_cross_yr);

    const auto t_relax =
        this->system_int_size() * t_cross_yr / (8 * std::log(this->system_int_size()));

    if (t_relax > 500'000)
        Logging::info("Relaxation Timescale at {} Myr", t_relax / 1'000'000);
    else
        Logging::info("Relaxation Timescale at {} yr", t_relax);

    Logging::info("==============================================================================");
}

auto System::particles_pos_at_step(uint idx)
    -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    std::vector<double> x, y, z;
    for (auto const &part : m_particles) {
        x.push_back(part.m_pos_history[idx].x());
        y.push_back(part.m_pos_history[idx].y());
        z.push_back(part.m_pos_history[idx].z());
    }
    return {x, y, z};
}

auto System::animate_particles() -> void {
    // gr.StartGIF("plots/particles.gif");
    const auto steps = m_particles[0].m_pos_history.size();
    mglGraph gr(0, 1440, 900);
    gr.SetRanges(-this->m_max_rad,
                 this->m_max_rad,
                 -this->m_max_rad,
                 this->m_max_rad,
                 -this->m_max_rad,
                 this->m_max_rad);

    for (uint i = 0; i < steps; i++) {
        Logging::info("Frame {}", i);
        gr.NewFrame(); // start frame
        gr.Rotate(50, 10);
        gr.Axis();

        const auto pos = particles_pos_at_step(i);
        const mglData x = std::get<0>(pos);
        const mglData y = std::get<1>(pos);
        const mglData z = std::get<2>(pos);

        gr.Dots(x, y, z, "r");
        gr.EndFrame(); // end frame
        std::stringstream ss;
        ss << "plots/animation/"
           << "frame-" << i << ".png";
        auto filename = ss.str();
        gr.WriteFrame(filename.c_str()); // save frame
        // gr.EndFrame();
    }
    // gr.CloseGIF();
}

auto System::get_analytic_mipd() const -> double {
    auto analytic_mipd = (((4 * std::numbers::pi) / 3) * this->m_half_mass_rad *
                          this->m_half_mass_rad * this->m_half_mass_rad);
    analytic_mipd /= this->system_int_size();
    analytic_mipd = std::pow(analytic_mipd, 1. / 3.);
    std::cerr << "DEBUGPRINT[1]: system.cpp:382: analytic_mipd=" << analytic_mipd << std::endl;
    return analytic_mipd;
}
