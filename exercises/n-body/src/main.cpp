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
        const auto rho_err = std_dev * System::k_non_dim_mass / shell.m_volume;

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
    gr.WriteJPEG("plots/hernquist.jpg");
    // gr.WritePNG("plots/hernquist.png");
    Logging::info("Hernquist plotted.");
}

auto plot_forces_step_2() {
    constexpr auto no_bins = 50;

// This should be a one time calculation and then save it in the System class
#if 0
    auto mean_inter_idst = g_system.precalc_mean_inter_part_dist();
    Logging::info("Mean inter particle distance: {}", mean_inter_idst);
#endif

    std::vector<double> analytic_force;
    std::vector<double> direct_force;
    std::vector<double> idx;

    // TODO: (dhub) Extract to System Class

    for (int i = 0; i <= no_bins; i++) {
        auto rad = g_system.convert_lin_to_log(no_bins, i);
        auto val = g_system.newton_force(rad);
        analytic_force.emplace_back(System::fit_log_to_plot(val));
    }

    // Initialize forces, with direct calculation
    g_system.calc_direct_initial_force();

    auto dir_force_hist = Histogram(no_bins, g_system, true);

    for (const auto &shell : dir_force_hist.m_shells) {

        // TODO: (aver) we need to convert vector force to the center of the spherical distribution
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
        // direct_force.emplace_back(val);
        idx.emplace_back(shell.m_lower_inc);

        // Logging::info("id: {}, {}", shell.m_lower_inc, val);
        // Logging::info(" {} particles in shell {}", shell.m_particles.size(), shell.m_index);
    }

    mglData x = idx;
    mglData aforce = analytic_force;
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
    gr.WriteJPEG("plots/forces.jpg");
    // gr.WritePNG("plots/forces.png");
}

auto plot_do_steps() {
    constexpr auto kTime = 10.;
    constexpr auto kDeltaTime = 0.001;
    // constexpr auto kFinalTime = kTime / kDeltaTime;

    for (double i = 0; i < kTime; i += kDeltaTime) {
        Logging::info("i = {}", i);
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
    // gr.WritePNG("plots/forces2.png");

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
}

auto tree_code() -> void {
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    double tolerance_angle = 0.5;
    TreeCode tree = TreeCode(root_cube, g_system.m_particles, tolerance_angle);
    tree.build();
    tree.tree_walk();

    mglGraph gr(0, 3000, 2000);

    // set plot parameters
    // gr.Rotate(50, 10); // Adjust for a better viewing angle
    // gr.SetRanges(-800, 800, -800, 800, -800, 800);
    // gr.Axis();
    // tree.plot(gr);
    // auto transform = g_system.transform_vectors();
    // mglData x = std::get<0>(transform);
    // mglData y = std::get<1>(transform);
    // mglData z = std::get<2>(transform);
    // gr.Dots(x,y,z, "r");
    // gr.WriteFrame("treecode.png");

    constexpr auto no_bins = 50;
    std::vector<double> numeric_force;
    std::vector<double> idx;

    g_system.m_particles = tree.m_particles;
    Logging::info("Creating Histogram with {} shells...", no_bins);
    auto dir_force_hist = Histogram(no_bins, g_system, true);

    for (const auto &shell : dir_force_hist.m_shells) {
        auto val = std::accumulate(
            shell.m_particles.begin(),
            shell.m_particles.end(),
            0.,
            [](double sum, const Particle3D &part) {
                auto norm = part.m_position.norm();
                auto projection = part.m_position.dot(part.m_tree_force) / norm;
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
    tree_code();

    Logging::info("Successfully quit!");
    return 0;
}
