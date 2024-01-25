#include "histogram.hpp"
#include "logging.hpp"
#include "node.hpp"
#include "particle.hpp"
#include "system.hpp"
#include "treecode.hpp"

#include "Eigen/Eigen"
#include "mgl2/mgl.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <format>
#include <numeric>
#include <string_view>

/// Particles read in and used in the tasks.
/// File Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

auto plot_rho_step_1() {
    std::vector<double> index;
    std::vector<double> hernquist_dens;
    std::vector<double> numeric_dens;
    std::vector<double> rho_error;

    constexpr auto no_bins = 50;
    const auto avg_parts = g_system.system_int_size() / no_bins;
    const auto std_dev = std::sqrt(avg_parts);

    // set minimal radius to something else than 0, otherwise errors ensue

    auto sum = 0.;
    Histogram hist(no_bins, g_system, true);

    for (const auto &shell : hist.m_shells) {

        auto hern_rho = g_system.density_hernquist((shell.m_lower_inc + shell.m_upper) / 2);

        // NOTE: (aver) \lambda = number of shells on average throughout all bins
        const auto rho_err = std_dev / shell.m_volume;
        // const auto rho_err = std_dev * shell.m_mass / shell.m_volume;

        // const auto k_no_parts_in_shell = shell.m_mass / Particle3D::km_non_dim_mass;
        // // NOTE: (aver) \lambda = number of shells in current bin
        // const auto rho_err =
        //     std::sqrt(k_no_parts_in_shell) * Particle3D::km_non_dim_mass / k_shell_volume;

        // // NOTE: (aver) \lambda = number of shells in a Hernquist bin
        // const auto no_parts_in_hern = (hern_rho * k_shell_volume) / Particle3D::km_non_dim_mass;
        // const auto rho_err =
        //     std::sqrt(no_parts_in_hern) * Particle3D::km_non_dim_mass / k_shell_volume;

        index.emplace_back(shell.m_lower_inc);

        hernquist_dens.emplace_back(System::fit_log_to_plot(hern_rho));
        numeric_dens.emplace_back(System::fit_log_to_plot(shell.m_density));
        rho_error.emplace_back(rho_err);

        sum += shell.m_mass;
    }

    assert(sum == 50010 || sum == 1001);

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
    gr.Axis();
    gr.SetCoor(mglLogLog);

    gr.Label('x', "Radius", 0);
    gr.Label('y', "Density", 0);

    gr.Plot(x, y_hern, "b");
    gr.AddLegend("Hernquist Density Profile", "b");

    gr.Plot(x, y_num, "r.");
    gr.AddLegend("Numeric Density Profile", "r.");

    gr.Error(x, y_num, y_err, "q");
    gr.AddLegend("Poissonian Error", "q");

    // finalize image
    gr.Legend();
    gr.WriteJPEG("plots/hernquist.jpg");
    gr.WritePNG("plots/hernquist.png");

    Logging::info("Hernquist plotted.");
}

auto plot_forces_step_2() {
    constexpr auto no_bins = 50;

// This should be a one time calculation and then save it in the System class
#if 0
    auto mean_inter_idst = g_system.precalc_mean_inter_part_dist();
    Logging::info("Mean inter-particle distance: {}", mean_inter_idst);

    auto analytic_mipd = (((4 * std::numbers::pi) / 3) * g_system.m_max_rad * g_system.m_max_rad *
                          g_system.m_max_rad);
    analytic_mipd /= g_system.system_int_size();
    analytic_mipd = std::pow(analytic_mipd, 1. / 3.);

    Logging::info("Analytic MIPD: {}", analytic_mipd);
#endif

    std::vector<double> analytic_force;

    // TODO: (dhub) Extract to System Class

    for (int i = 0; i <= no_bins; i++) {
        auto rad = g_system.convert_lin_to_log(no_bins, i);
        auto val = g_system.newton_force(rad);
        analytic_force.emplace_back(System::fit_log_to_plot(val));
    }

    mglData aforce = analytic_force;

    for (const auto &div : System::softening_divisors) {
        std::vector<double> direct_force;
        std::vector<double> idx;

        System::s_softening = System::s_softening_length / div;
        // Initialize forces, with direct calculation
        g_system.calc_direct_initial_force();

        auto dir_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : dir_force_hist.m_shells) {

            // TODO: (aver) we need to convert vector force to the center of the spherical
            // distribution
            auto val = std::accumulate(shell.m_particles.begin(),
                                       shell.m_particles.end(),
                                       0.,
                                       [](double sum, const Particle3D &part) {
                                           auto norm = part.m_position.norm();
                                           auto projection =
                                               part.m_position.dot(part.m_direct_force) / norm;
                                           return sum + projection;
                                       });
            val = shell.shell_int_size() == 0 ? 0. : val / shell.shell_int_size();

            direct_force.emplace_back(System::fit_log_to_plot(val));
            idx.emplace_back(shell.m_lower_inc);

        }

        mglData x = idx;
        mglData nforce = direct_force;

        auto y_min = std::min(aforce.Minimal(), nforce.Minimal());
        auto y_max = std::max(aforce.Maximal(), nforce.Maximal());

        mglGraph gr(0, 3000, 2000);

        gr.SetRange('x', x);
        gr.SetRange('y', y_min, y_max);
        // gr.SetRange('y', nforce);

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
        gr.WriteJPEG(std::format("plots/forces_{}.jpg", div).c_str());
        // gr.WritePNG(std::format("plots/forces_{}.png", div).c_str());
    }
}

auto plot_do_steps() {
    System::s_softening = System::s_softening_length / 200;
    g_system.calc_direct_initial_force();

    const auto kTime = g_system.m_t_cross * 5;
    constexpr auto kDeltaTime = 0.0000001;
    // constexpr auto kFinalTime = kTime / kDeltaTime;

    for (double t = 0; t < kTime; t += kDeltaTime) {
        Logging::info("t = {}", t);
        g_system.solver_do_step(kDeltaTime);
    }

    std::vector<double> numeric_dens;
    std::vector<double> direct_force;
    std::vector<double> idx;

    constexpr auto no_bins = 50;
    Logging::info("Creating Histogram with {} shells...", no_bins);
    auto dir_force_hist = Histogram(no_bins, g_system, true);

    for (const auto &shell : dir_force_hist.m_shells) {
        auto val = std::accumulate(shell.m_particles.begin(),
                                   shell.m_particles.end(),
                                   0.,
                                   [](double sum, const Particle3D &part) {
                                       auto norm = part.m_position.norm();
                                       auto projection =
                                           part.m_position.dot(part.m_direct_force) / norm;
                                       return sum + projection;
                                   });
        val = shell.shell_int_size() == 0 ? 0. : val / shell.shell_int_size();

        numeric_dens.emplace_back(System::fit_log_to_plot(shell.m_density));
        direct_force.emplace_back(System::fit_log_to_plot(val));
        idx.emplace_back(shell.m_lower_inc);
    }

    mglData x = idx;
    mglData nforce = direct_force;
    mglData ndens = numeric_dens;

    mglGraph gr(0, 3000, 2000);

    gr.SetRange('x', x);
    gr.SetRange('y', nforce);

    gr.SetFontSize(2);
    gr.SetCoor(mglLogX);
    gr.Axis();

    gr.Label('x', "Radius [l]", 0);
    gr.Label('y', "Force", 0);

    gr.Plot(x, nforce, "r.");
    gr.AddLegend("Numeric", "r.");

    gr.Legend();
    gr.WriteJPEG("plots/forces2.jpg");
    gr.WritePNG("plots/forces2.png");

    gr.ClearFrame();
    gr.ClearLegend();
    gr.ResetFrames();

    gr.SetRange('x', x);
    gr.SetRange('y', ndens);

    gr.SetFontSize(2);
    gr.SetCoor(mglLogLog);
    gr.Axis();

    gr.Plot(x, ndens, "r.");
    gr.AddLegend("Numeric", "r.");

    gr.Legend();
    gr.WriteJPEG("plots/density2.jpg");
    gr.WritePNG("plots/density2.png");
}

auto tree_code() -> void {
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    double tolerance_angle = 0.5;
    TreeCode tree = TreeCode(root_cube, g_system.m_particles, tolerance_angle);
    tree.build();
    // tree.tree_walk();

    mglGraph gr(0, 3000, 2000);

    // set plot parameters
    gr.Rotate(50, 10); // Adjust for a better viewing angle
    gr.SetRanges(-800, 800, -800, 800, -800, 800);
    gr.Axis();
    tree.plot(gr);
    auto transform = g_system.transform_vectors();
    mglData x = std::get<0>(transform);
    mglData y = std::get<1>(transform);
    mglData z = std::get<2>(transform);
    gr.Dots(x, y, z, "r");
    gr.WriteFrame("plots/treecode.png");

    /*
    constexpr auto no_bins = 50;
    std::vector<double> numeric_force;
    std::vector<double> idx;

    g_system.m_particles = tree.m_particles;
    Logging::info("Creating Histogram with {} shells...", no_bins);
    auto dir_force_hist = Histogram(no_bins, g_system, true);

    for (const auto &shell : dir_force_hist.m_shells) {
        auto val = std::accumulate(shell.m_particles.begin(),
                                   shell.m_particles.end(),
                                   0.,
                                   [](double sum, const Particle3D &part) {
                                       auto norm = part.m_position.norm();
                                       auto projection =
                                           part.m_position.dot(part.m_tree_force) / norm;
                                       return sum + projection;
                                   });
        val = shell.shell_int_size() == 0 ? 0. : val / shell.shell_int_size();

        numeric_force.emplace_back(System::fit_log_to_plot(val));
        idx.emplace_back(shell.m_lower_inc);
    }

    mglData x = idx;
    mglData nforce = numeric_force;

    gr.SetRange('x', x);
    gr.SetRange('y', nforce);

    gr.SetFontSize(2);
    gr.SetCoor(mglLogX);
    gr.Axis();

    gr.Label('x', "Radius [l]", 0);
    gr.Label('y', "Force", 0);

    gr.Plot(x, nforce, "r.");
    gr.AddLegend("Numeric", "r.");

    gr.Legend();
    gr.WriteJPEG("plots/forces3.jpg");
    gr.WritePNG("plots/forces3.png");
    */
}

/// Calculate timescales with remultiplied factors and units
auto calc_real_relaxation() {
    constexpr auto parsec_to_km_factor = 3.0857e+13;
    Logging::info("==============================================================================");
    Logging::info("Calculating relaxation timescale");

    auto rel_hist = Histogram(100'000, g_system);
    for (auto &shell : rel_hist.m_shells) {
        shell.m_mass = shell.shell_int_size() * System::k_dim_mass;
        shell.update_density();
    }

    const auto r_hm = g_system.calc_half_mass_radius(rel_hist.m_shells);
    const auto M = g_system.system_int_size() * System::k_dim_mass;
    Logging::info("Half Mass Radius at     {} pc", r_hm);

    // WARN: (aver) What to do? Use G=1 or use the real unit of G?
    const auto v_c = std::sqrt(System::k_G * M * 0.5 / r_hm);
    Logging::info("Circular Velocity at    {} km / s", v_c);

    // const auto t_cross = r_hm / v_c;
    // Logging::info("Crossing Timescale at   {} pc / (km/s)", t_cross);
    const auto t_cross_km = (r_hm * parsec_to_km_factor) / v_c;
    const auto t_cross_yr = t_cross_km / (60. * 60. * 24. * 365.25);
    Logging::info("Crossing Timescale at   {} yr", t_cross_yr);

    const auto t_relax =
        g_system.system_int_size() / (8 * std::log(g_system.system_int_size()) * t_cross_yr);
    Logging::info("Relaxation Timescale at {} yr", t_relax);

    Logging::info("==============================================================================");
}

auto main(const int argc, const char *const argv[]) -> int {
    if (argc != 2) {
        Logging::err("Please supply a single file argument!");
        return -1;
    }
    // initialize the g_system variable
    g_system.init_system(argv[1]);

    // ============================================================================================
    // task 1
    // ============================================================================================
    plot_rho_step_1();
    // TODO: (aver)
    // - We still need to do a comparison between different softeining values and discuss their
    //      significance
    // - Also explain dependence of force calculation on direct force calculation
    plot_forces_step_2();

    // ============================================================================================
    // task 2
    // ============================================================================================
    // plot_do_steps();
    // tree_code();

    calc_real_relaxation();
    Logging::info("Successfully quit!");
    return 0;
}
