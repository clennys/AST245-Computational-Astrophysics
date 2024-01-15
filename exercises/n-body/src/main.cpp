#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "system.hpp"

#include <cassert>
#include <cmath>
#include <format>
#include <mgl2/mgl.h>
#include <string_view>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

// TODO: (aver) consider making this a factory like static function on the `System` class
auto init_system(const std::string_view &path) {
    g_system = System(path);
    g_system.precalc_consts();

    auto furthest_particle = g_system.get_max_distance();

    // use this shell to calculcate half_mass_radius, and scale_length
    Histogram hist(100'000, furthest_particle.m_distance, g_system);

    // g_system.update_half_mass(hist.m_shells);
    g_system.update_half_mass_radius(hist.m_shells);
    g_system.update_scale_length();
    g_system.m_softening = g_system.m_max_rad / std::pow(g_system.m_total_mass, 1. / 3.);
    Particle3D::s_softening = g_system.m_softening;

    Logging::info("Total mass of system:       {:<12}", g_system.m_total_mass);
    Logging::info("Half mass radius of system: {:>12.10f}", g_system.m_half_mass_rad);
    Logging::info("Scaling length of system:   {:>12.10f}", g_system.m_scale_length);
    Logging::info("Softening of system:        {:>12.10f}", g_system.m_softening);
}

auto plot_rho_step_1() {
    std::vector<double> index;
    std::vector<double> hernquist_dens;
    std::vector<double> numeric_dens;
    std::vector<double> rho_error;

    // precalculated (by the compiler!) constant for shell volume
    constexpr auto k_shell_vol_pref = 4. / 3. * std::numbers::pi;
    constexpr auto no_bins = 50;
    const auto avg_parts = static_cast<int>(g_system.m_particles.size()) / no_bins;
    const auto std_dev = std::sqrt(avg_parts);

    // set minimal radius to something else than 0, otherwise errors ensue

    auto sum = 0.;
#if 1
    auto furthest_particle = g_system.get_max_distance();
    Histogram hist(no_bins, furthest_particle.m_distance, g_system, true);
    for (const auto &shell : hist.m_shells) {
        // FIXME: (aver) fix plotting of hernquist

        Logging::info(
            "bin: {}, lower: {}, upper {}", shell.m_index, shell.m_lower_inc, shell.m_upper);
        Logging::info("Shell mass: {}", shell.m_mass);

        // TODO: (aver) move this calculation into a method of shell
        const auto k_shell_volume =
            k_shell_vol_pref * (std::pow(shell.m_upper, 3) - std::pow(shell.m_lower_inc, 3));

        // TODO: (aver) move this into a method of shell as well
        const auto k_shell_rho = shell.m_mass / k_shell_volume;
        std::cerr << k_shell_rho << std::endl;

        auto hern_rho = g_system.density_hernquist((shell.m_lower_inc + shell.m_upper) / 2);
        std::cerr << hern_rho << std::endl;

        // NOTE: (aver) \lambda = number of shells on average throughout all bins
        const auto rho_err = std_dev * Particle3D::km_non_dim_mass / k_shell_volume;

        // const auto k_no_parts_in_shell = shell.m_mass / Particle3D::km_non_dim_mass;
        // // NOTE: (aver) \lambda = number of shells in current bin
        // const auto rho_err =
        //     std::sqrt(k_no_parts_in_shell) * Particle3D::km_non_dim_mass / k_shell_volume;

        // // NOTE: (aver) \lambda = number of shells in a Hernquist bin
        // const auto no_parts_in_hern = (hern_rho * k_shell_volume) / Particle3D::km_non_dim_mass;
        // const auto rho_err =
        //     std::sqrt(no_parts_in_hern) * Particle3D::km_non_dim_mass / k_shell_volume;

        // index.emplace_back(static_cast<int>(shell.m_index));
        index.emplace_back(shell.m_lower_inc);

        hernquist_dens.emplace_back(System::fit_log_to_plot(hern_rho));
        numeric_dens.emplace_back(System::fit_log_to_plot(k_shell_rho));
        rho_error.emplace_back(rho_err);

        sum += shell.m_mass;
        // Logging::info("\n\tH: {}\n\tN: {}", hern_rho, k_shell_rho);
    }

#else
    for (int r = 0; r <= no_bins; r++) {
        auto lower_rad = g_system.convert_lin_to_log(no_bins, r);
        auto upper_rad = g_system.convert_lin_to_log(no_bins, r + 1);
        Logging::info("bin: {}, lower: {}, upper {}", r, lower_rad, upper_rad);

        const auto k_shell_mass = g_system.get_constrained_shell_mass(lower_rad, upper_rad);
        Logging::info("Shell mass: {}", k_shell_mass);
        sum += k_shell_mass;
        const auto k_shell_volume =
            k_shell_vol_pref * (std::pow(upper_rad, 3) - std::pow(lower_rad, 3));

        const auto k_no_parts_in_shell = k_shell_mass / Particle3D::km_non_dim_mass;

        const auto k_shell_rho = k_shell_mass / k_shell_volume;
        const auto k_hern_rho = g_system.density_hernquist((lower_rad + upper_rad) / 2);

        // NOTE: (aver) \lambda = number of shells in current bin
        // const auto k_rho_error =
        //     std::sqrt(k_no_parts_in_shell) * Particle3D::km_non_dim_mass / k_shell_volume;

        // NOTE: (aver) \lambda = number of shells in a Hernquist bin
        const auto no_parts_in_hern = (k_hern_rho * k_shell_volume) / Particle3D::km_non_dim_mass;
        // const auto k_rho_error =
        //     std::sqrt(no_parts_in_hern) * Particle3D::km_non_dim_mass / k_shell_volume;

        // NOTE: (aver) \lambda = number of shells on average throughout all bins
        const auto k_rho_error = std_dev * Particle3D::km_non_dim_mass / k_shell_volume;

        index.emplace_back(lower_rad);
        hernquist_dens.emplace_back(System::fit_log_to_plot(k_hern_rho));
        numeric_dens.emplace_back(System::fit_log_to_plot(k_shell_rho));
        rho_error.emplace_back(k_rho_error);

        sum += shell.m_mass;
    }
#endif

    // Logging::dbg("total mass {}, expected 50010", sum);
    assert(sum == 50010);

    mglData x = index;
    mglData y_hern = hernquist_dens;
    mglData y_num = numeric_dens;
    mglData y_err = rho_error;

    auto y_min = std::min(y_hern.Minimal(), y_num.Minimal());
    auto y_max = std::max(y_hern.Maximal(), y_num.Maximal());

    // set plot parameters
    mglGraph gr(0, 3000, 2000);

    gr.SetRange('x', x);
    gr.SetRange('y', y_min, y_max);

    gr.SetFontSize(2);
    gr.SetCoor(mglLogLog);
    gr.Axis();

    gr.Label('x', "Radius [l]", 0);
    gr.Label('y', "Density [m]/[l]^3", 0);

    gr.Plot(x, y_hern, "b");
    gr.AddLegend("Hernquist Density Profile", "b");

    gr.Plot(x, y_num, "r .");
    gr.AddLegend("Numeric Density Profile", "r .");

    gr.Error(x, y_num, y_err, "q");
    gr.AddLegend("Poissonian Error", "q");

    gr.Legend();
    gr.WriteFrame("hernquist.jpg");
    gr.WriteFrame("hernquist.png");
    Logging::info("Hernquist plotted.");
}

auto plot_forces_step_2() {

    constexpr auto no_bins = 100;

    std::vector<double> analytic_force;
    std::vector<double> direct_force;
    std::vector<double> idx;

    // TODO: (dhub) Extract to System Class

    for (int i = 0; i <= no_bins; i++) {
        auto val = g_system.newton_force(g_system.convert_lin_to_log(no_bins, i));
        Logging::info("{}", val);
        analytic_force.emplace_back(System::fit_log_to_plot(val));
    }

    auto furthest_particle = g_system.get_max_distance();
    g_system.calc_direct_force();
    auto dir_force_hist = Histogram(no_bins, furthest_particle.m_distance, g_system, true);
    for (const auto &shell : dir_force_hist.m_shells) {
        // TODO: (aver) add logic for getting values from direct force calculation

        auto avg_force = Eigen::Vector3d({0., 0., 0.});
        for (const auto &part : shell.m_particles) {
            avg_force += part.m_direct_force;
        }

        auto val = avg_force.norm() == 0.
                       ? 0.
                       : avg_force.norm() / static_cast<int>(shell.m_particles.size());
        direct_force.emplace_back(val);
        // Logging::info("id: {}, {}", shell.m_index, val);

        idx.emplace_back(shell.m_lower_inc);
    }

    mglData x = idx;
    mglData aforce = analytic_force;
    mglData nforce = direct_force;

    auto y_min = std::min(aforce.Minimal(), nforce.Minimal());
    auto y_max = std::max(aforce.Maximal(), nforce.Maximal());

    mglGraph gr(0, 3000, 2000);

    gr.SetRange('x', x);
    // gr.SetRange('y', aforce);
    gr.SetRange('y', y_min, y_max);

    gr.SetFontSize(2);
    gr.SetCoor(mglLogX);
    gr.Axis();

    gr.Label('x', "Radius [l]", 0);
    gr.Label('y', "Force", 0);

    gr.Plot(x, aforce, "b");
    gr.AddLegend("Analytic", "b");

    gr.Plot(x, nforce, "r.");
    gr.AddLegend("Numeric", "r.");

    gr.Legend();
    gr.WritePNG("forces.png");
    gr.WritePNG("forces.jpg");
}

auto main(const int argc, const char *const argv[]) -> int {
    if (argc != 2) {
        Logging::err("Please supply a single file argument!");
        return -1;
    }
    // initialize the g_system variable
    init_system(argv[1]);

    // task 1
    plot_rho_step_1();
    // plot_forces_step_2();

    // task 2

    Logging::info("Successfully quit!");
    return 0;
}
