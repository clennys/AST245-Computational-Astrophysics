#include "histogram.hpp"
#include "logging.hpp"
#include "node.hpp"
#include "particle.hpp"
#include "system.hpp"
#include "treecode.hpp"

#include "Eigen/Eigen"
#include "mgl2/mgl.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <string_view>

/// Particles read in and used in the tasks.
/// File Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

/// Plots the density distribution numerically compared with the analytical distribution from the
/// Hernquist paper
auto plot_rho_step_1() {
    std::vector<double> index;
    std::vector<double> hernquist_dens;
    std::vector<double> numeric_dens;
    std::vector<double> rho_error;

    constexpr auto no_bins = 50;
    const auto avg_parts = g_system.system_int_size() / no_bins;
    const auto std_dev = std::sqrt(avg_parts);

    // set minimal radius to something else than 0, otherwise errors ensue

    auto hist = Histogram(no_bins, g_system, true);

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
    }

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

/// Plots the forces calculated by direct summation compared to the analytical force distribution
/// from the Hernquist paper
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
        g_system.precalc_direct_initial_force();

        auto dir_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : dir_force_hist.m_shells) {
            auto val = shell.get_avg_direct_force();

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

/// Plot a final step after a predefined time integration
auto plot_do_steps() {
    System::s_softening = System::s_softening_length / 200;
    g_system.precalc_direct_initial_force();

    constexpr auto tau = 0.01;
    const auto kTime = g_system.m_t_cross * 3;
    const auto kDeltaTime = tau * g_system.m_t_cross;

    for (double t = 0.; t < kTime; t += kDeltaTime) {
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
        auto val = shell.get_avg_direct_force();

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

/// Create an octree and calculate the forces by running multipole expansion on the tree
auto tree_code() -> void {
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    double tolerance_angle = 0.1;
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
        auto val = shell.get_avg_tree_force();

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

/// Creates a GIF showing the variances of the force distribution integrated via tree
auto plot_gif_steps_tree() {
    constexpr auto no_bins = 50;
    constexpr auto tau = 0.01;
    const auto kTime = g_system.m_t_cross * 2;
    const auto kDeltaTime = tau * g_system.m_t_cross;

    g_system.reset_system();
    System::s_softening = System::s_softening_length / 200;
    g_system.precalc_direct_initial_force();

    constexpr auto tolerance_angle = 0.1;
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    auto tree = TreeCode(root_cube, g_system.m_particles, tolerance_angle);

    mglGraph gr(0, 1440, 900);
    gr.StartGIF("steps_tree_forces.gif");

    tree.build();
    tree.tree_walk();

    for (double t = 0.; t < kTime; t += kDeltaTime) {
        Logging::info("t = {}", t);
        gr.NewFrame();

        // std::vector<double> numeric_dens;
        std::vector<double> tree_force;
        std::vector<double> idx;

        g_system.m_particles = tree.m_particles;
        auto dir_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : dir_force_hist.m_shells) {
            auto val = shell.get_avg_tree_force();

            // numeric_dens.emplace_back(System::fit_log_to_plot(shell.m_density));
            tree_force.emplace_back(System::fit_log_to_plot(val));
            idx.emplace_back(shell.m_lower_inc);
        }

        mglData x = idx;
        mglData nforce = tree_force;
        // mglData ndens = numeric_dens;

        gr.SetRange('x', x);
        gr.SetRange('y', nforce);

        gr.SetFontSize(2);
        gr.SetCoor(mglLogX);
        gr.Axis();

        gr.Label('x', "Radius [l]", 0);
        gr.Label('y', "Force", 0);

        gr.Plot(x, nforce, "r.");
        gr.AddLegend("Numeric", "r.");

        // Add a text box at the top right
        gr.Title(std::format("t: {}", t).c_str());

        gr.Legend();

        gr.EndFrame();

        tree.reset_tree();
        tree.build();
        tree.tree_step(kDeltaTime);
    }

    gr.CloseGIF();
}

/// Creates a GIF showing the variances of the force distribution integrated via direct calculation
auto plot_gif_steps() {
    constexpr auto no_bins = 50;
    constexpr auto tau = 0.01;
    const auto kTime = g_system.m_t_cross * 2;
    const auto kDeltaTime = tau * g_system.m_t_cross;

    g_system.reset_system();
    System::s_softening = System::s_softening_length / 200;
    g_system.precalc_direct_initial_force();

    mglGraph gr(0, 1440, 900);
    gr.StartGIF("steps_forces.gif");

    for (double t = 0.; t < kTime; t += kDeltaTime) {
        Logging::info("t = {}", t);
        gr.NewFrame();

        std::vector<double> numeric_dens;
        std::vector<double> direct_force;
        std::vector<double> idx;

        auto dir_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : dir_force_hist.m_shells) {
            auto val = shell.get_avg_direct_force();

            numeric_dens.emplace_back(System::fit_log_to_plot(shell.m_density));
            direct_force.emplace_back(System::fit_log_to_plot(val));
            idx.emplace_back(shell.m_lower_inc);
        }

        Logging::info("t = {}", t);
        g_system.solver_do_step(kDeltaTime);

        mglData x = idx;
        mglData nforce = direct_force;
        mglData ndens = numeric_dens;

        gr.SetRange('x', x);
        gr.SetRange('y', nforce);

        gr.SetFontSize(2);
        gr.SetCoor(mglLogX);
        gr.Axis();

        gr.Label('x', "Radius [l]", 0);
        gr.Label('y', "Force", 0);

        gr.Plot(x, nforce, "r.");
        gr.AddLegend("Numeric", "r.");

        // Add a text box at the top right
        gr.Title(std::format("t: {}", t).c_str());

        gr.Legend();

        gr.EndFrame();
    }

    gr.CloseGIF();
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
    // plot_rho_step_1();
    // TODO: (aver)
    // - We still need to do a comparison between different softeining values and discuss their
    //      significance
    // - Also explain dependence of force calculation on direct force calculation
    // plot_forces_step_2();

    // g_system.calc_real_relaxation();
    // ============================================================================================
    // task 2
    // ============================================================================================

    // plot_do_steps();
    // tree_code();
    plot_gif_steps();

    Logging::info("Successfully quit!");
    return 0;
}
