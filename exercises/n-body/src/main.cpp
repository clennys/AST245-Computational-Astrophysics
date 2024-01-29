#include "histogram.hpp"
#include "logging.hpp"
#include "node.hpp"
#include "particle.hpp"
#include "system.hpp"
#include "timer.hpp"
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
    auto hist = Histogram(no_bins, g_system, true);

    for (const auto &shell : hist.m_shells) {

        auto hern_rho = g_system.density_hernquist((shell.m_lower_inc + shell.m_upper) / 2.);

        // NOTE: (aver) \lambda = number of shells on average throughout all bins
        const auto rho_err = std_dev * System::k_dim_mass / shell.m_volume;

        // // NOTE: (aver) \lambda = number of shells in current bin
        // const auto rho_err =
        //     std::sqrt(shell.shell_int_size()) * System::k_dim_mass / shell.m_volume;

        // // NOTE: (aver) \lambda = number of shells in a Hernquist bin
        // const auto no_parts_in_hern = (hern_rho * shell.m_volume) / System::k_dim_mass;
        // const auto rho_err =
        //     std::sqrt(no_parts_in_hern) * System::k_dim_mass / shell.m_volume;

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
    gr.SetFontSize(2);

    gr.SetRange('x', x);
    gr.SetRange('y', y_min, y_max);
    gr.SetCoor(mglLogLog);
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
    gr.WritePNG("plots/png/hernquist.png");

    Logging::info("Hernquist plotted.");
}

/// Plots the forces calculated by direct summation compared to the analytical force distribution
/// from the Hernquist paper
auto plot_forces_step_2() {
    constexpr auto no_bins = 50;

// This should be a one time calculation and then save it in the System class
#if 1
    auto mean_inter_idst = g_system.precalc_mean_inter_part_dist();
    Logging::info("Mean inter-particle distance: {}", mean_inter_idst);
    // use the half mass radius as requested in the instructions
    System::s_softening_length = g_system.get_analytic_mipd();
    Logging::info("Analytic MIPD: {}", g_system.get_analytic_mipd());
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

        // System::s_softening = System::s_softening_length / div;
        System::s_softening = g_system.get_analytic_mipd() / div;

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

        gr.Title(std::format("Softening: {}", System::s_softening ).c_str());
        gr.Legend();
        gr.WriteJPEG(std::format("plots/forces_{}.jpg", div).c_str());
        gr.WritePNG(std::format("plots/png/forces_{}.png", div).c_str());
    }
}

/// Create an octree and calculate the forces by running multipole expansion on the tree
auto tree_code() -> void {
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    auto tolerance_angles = {0.1, 0.3, 0.5, 0.7, 0.9};
    for (auto tol : tolerance_angles) {
        TreeCode tree = TreeCode(root_cube, g_system.m_particles, tol);
        tree.build();
        // auto cost = tree.computational_cost(g_system.m_half_mass_rad,
        // g_system.k_mean_inter_dist); std::cerr << "DEBUGPRINT[2]: main.cpp:251: cost=" << cost <<
        // std::endl;

        mglGraph gr_tree(0, 3000, 2000);
        gr_tree.SetFontSize(2);

        // set plot parameters
        gr_tree.Rotate(50, 10); // Adjust for a better viewing angle
        gr_tree.SetRanges(-800, 800, -800, 800, -800, 800);
        gr_tree.Axis();
        tree.plot(gr_tree);
        auto transform = g_system.transform_vectors();
        mglData x = std::get<0>(transform);
        mglData y = std::get<1>(transform);
        mglData z = std::get<2>(transform);
        gr_tree.Dots(x, y, z, "r");
        gr_tree.WriteFrame(std::format("plots/treecode_{}.jpg", tol).c_str());
        gr_tree.WriteFrame(std::format("plots/png/treecode_{}.png", tol).c_str());

        constexpr auto no_bins = 50;
        std::vector<double> tree_force;
        std::vector<double> direct_force;
        std::vector<double> analytic_force;
        std::vector<double> idx;

        for (int i = 0; i <= no_bins; i++) {
            auto rad = g_system.convert_lin_to_log(no_bins, i);
            auto val = g_system.newton_force(rad);
            analytic_force.emplace_back(System::fit_log_to_plot(val));
        }

        auto analytic_mipd = (((4 * std::numbers::pi) / 3) * g_system.m_half_mass_rad *
                              g_system.m_half_mass_rad * g_system.m_half_mass_rad);
        analytic_mipd /= g_system.system_int_size() * .5;
        analytic_mipd = std::pow(analytic_mipd, 1. / 3.);

        System::s_softening = analytic_mipd / 200;
        g_system.precalc_direct_initial_force();

        auto dir_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : dir_force_hist.m_shells) {
            auto df = shell.get_avg_direct_force();
            direct_force.emplace_back(System::fit_log_to_plot(df));
        }

        tree.tree_walk();
        auto err = tree.m_force_error;
        std::cerr << "DEBUGPRINT[1]: main.cpp:253: err=" << err << std::endl;
        g_system.m_particles = tree.m_particles;
        auto tree_force_hist = Histogram(no_bins, g_system, true);

        for (const auto &shell : tree_force_hist.m_shells) {
            auto tf = shell.get_avg_tree_force();
            tree_force.emplace_back(System::fit_log_to_plot(tf));
            idx.emplace_back(shell.m_lower_inc);
        }

        mglGraph gr(0, 3000, 2000);
        gr.SetFontSize(2);

        mglData x_2 = idx;
        mglData tree_y = tree_force;
        mglData direct_y = direct_force;
        mglData hern_y = analytic_force;

        const auto y_min = std::min(tree_y.Minimal(), direct_y.Minimal());
        const auto y_max = std::max(tree_y.Maximal(), direct_y.Maximal());

        gr.SetRange('x', x_2);
        gr.SetRange('y', y_min, y_max);

        gr.SetCoor(mglLogX);
        gr.Axis();

        gr.Label('x', "Radius", 0);
        gr.Label('y', "Force", 0);

        gr.Plot(x_2, tree_y, "r.");
        gr.AddLegend(std::format("Tree Force with \\theta: {}", tol).c_str(), "r.");

        gr.Plot(x_2, direct_y, "b.");
        gr.AddLegend(std::format("Direct Force with \\epsilon: {:e}", System::s_softening).c_str(),
                     "b.");

        gr.Plot(x_2, hern_y, "q");
        gr.AddLegend(std::format("Hernquist Analytical Force", System::s_softening).c_str(), "q");

        gr.Legend();
        gr.WriteFrame(std::format("plots/tree_force_{}.jpg", tol).c_str());
        gr.WriteFrame(std::format("plots/png/tree_force_{}.png", tol).c_str());

        tree.reset_tree();
    }
}

/// Creates a GIF showing the variances of the force distribution integrated via tree
auto plot_gif_steps_tree() {
    constexpr static auto no_bins = 50;
    constexpr static auto t_factor = 2;
    constexpr auto eta = 0.01;
    const auto kTime = g_system.m_t_cross * t_factor;
    const auto kDeltaTime = eta * g_system.m_t_cross;
    const auto no_steps = kTime / kDeltaTime;
    std::vector<double> en_tot;

    // g_system.reset_system();
    // System::s_softening = System::s_softening_length / 200;
    // g_system.precalc_direct_initial_force();

    constexpr auto tolerance_angle = 0.1;
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    auto tree = TreeCode(root_cube, g_system.m_particles, tolerance_angle);

    mglGraph gr(0, 1440, 900);
    gr.StartGIF(std::format("tree_forces_{}_eta_{}_soft.gif", t_factor, eta).c_str());

    tree.build();
    tree.tree_walk();

    int step = 0;
    for (double t = 0.; t < kTime; t += kDeltaTime) {
        Logging::info("step: {}/{}, t: {:e}", step, no_steps, t);
        gr.NewFrame();
        tree.total_energy(System::k_dim_mass);
        auto pot = tree.m_pot_energy;
        std::cerr << "DEBUGPRINT[10]: main.cpp:297: pot=" << pot << std::endl;
        auto kin = tree.m_kin_energy;
        std::cerr << "DEBUGPRINT[11]: main.cpp:299: kin=" << kin << std::endl;
        auto tot = tree.m_tot_energy;
        en_tot.push_back(tot);
        std::cerr << "DEBUGPRINT[12]: main.cpp:301: tot=" << tot << std::endl;

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

        gr.Label('x', "Radius", 0);
        gr.Label('y', "Force", 0);

        gr.Plot(x, nforce, "r.");
        gr.AddLegend(std::format("Tree Force with \\theta: {}", tolerance_angle).c_str(), "r.");

        // Add a text box at the top right
        gr.Title(std::format("Frame: {}\nt: {:e}", step, t).c_str());
        gr.Legend();
        gr.EndFrame();

        tree.reset_tree();
        tree.build();
        tree.tree_step(kDeltaTime);
        step++;
    }

    gr.CloseGIF();
    return en_tot;
}

/// Creates a GIF showing the variances of the force distribution integrated via direct calculation
auto plot_gif_steps(double eta, double div)
    -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    constexpr static auto no_bins = 50;
    constexpr static auto t_factor = 5;
    const auto kTime = g_system.m_t_cross * t_factor;
    const auto kDeltaTime = eta * g_system.m_t_cross;
    const auto no_steps = kTime / kDeltaTime;
    std::vector<double> tot_en, t_steps, virial;

    auto analytic_mipd = (((4 * std::numbers::pi) / 3) * g_system.m_half_mass_rad *
                          g_system.m_half_mass_rad * g_system.m_half_mass_rad);
    analytic_mipd /= g_system.system_int_size() * .5;
    analytic_mipd = std::pow(analytic_mipd, 1. / 3.);

    System::s_softening_length = analytic_mipd;

    g_system.reset_system();
    System::s_softening = g_system.get_analytic_mipd() / div;
    g_system.precalc_direct_initial_force();

    mglGraph gr(0, 1440, 900);
    gr.StartGIF(std::format("steps_forces_{}_eta_{}_soft_{}.gif", t_factor, eta, div).c_str());
    gr.SetFontSize(2);

    int step = 1;
    for (double t = 0.; t < kTime; t += kDeltaTime) {
        Logging::info("step: {}/{}, t: {:e}", step, no_steps, t);
        gr.NewFrame();
        g_system.total_energy();
        auto tot = g_system.m_tot_energy;
        auto vir = g_system.m_delta_energy;
        tot_en.push_back(tot);
        t_steps.push_back(t);
        virial.push_back(vir);
        std::cerr << "DEBUGPRINT[12]: main.cpp:301: tot=" << tot << std::endl;
        std::cerr << "DEBUGPRINT[1]: main.cpp:388: delta e=" << vir << std::endl;

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

        g_system.solver_do_step(kDeltaTime);

        mglData x = idx;
        mglData nforce = direct_force;
        mglData ndens = numeric_dens;

        gr.SetRange('x', x);
        gr.SetRange('y', nforce);

        gr.SetCoor(mglLogX);
        gr.Axis();

        gr.Label('x', "Radius", 0);
        gr.Label('y', "Force", 0);

        gr.Plot(x, nforce, "r.");
        gr.AddLegend(std::format("Direct Force with Softening: {:e}", System::s_softening).c_str(),
                     "r.");

        // Add a text box at the top right
        gr.Title(std::format("Frame: {}\nt: {:e}", step, t).c_str());

        gr.Legend();

        gr.EndFrame();
        step++;
    }
    gr.CloseGIF();
    return {t_steps, tot_en, virial};
}

auto plot_en(const std::vector<double> &t_step,
             const std::vector<double> &en_direct,
             double eta,
             double div,
             bool vir = false) {

    mglGraph gr(0, 1440, 900);
    mglData x = t_step;
    mglData y1 = en_direct;

    // double max_value = *std::max_element(en_direct.begin(), en_direct.end());

    gr.SetRange('x', x);
    // gr.SetRange('y', 0, max_value);
    gr.SetRange('y', y1);

    gr.SetFontSize(2);
    gr.Axis();

    gr.Label('x', "Time Step", 0);
    if (vir) {
        gr.Label('y', "\\Delta E", 0);
    } else {
        gr.Label('y', "E_{tot}", 0);
    }

    gr.Plot(x, y1, "r.");
    gr.AddLegend("Direct Sum", "r.");

    // Add a text box at the top right
    gr.Legend();
    if (vir) {
        gr.Title(std::format("Change in Total Energy with \\eta={}, div={}", eta, div).c_str());
        gr.WritePNG(std::format("plots/png/change_energy_{}_{}.png", eta, div).c_str());

    } else {
        gr.Title(std::format("Numerical Relaxation with eta={}, div={}", eta, div).c_str());
        gr.WritePNG(std::format("plots/png/num_relax_en_{}_{}.png", eta, div).c_str());
    }
}
auto comp_cost() -> void {
    BoundingCube root_cube = g_system.calc_overall_bounding_cube();
    double tolerance_angle = 1.0;
    TreeCode tree = TreeCode(root_cube, g_system.m_particles, tolerance_angle);
    tree.build();
    tree.tree_walk();
    auto mp = tree.m_count_mp_inter;
    std::cerr << "DEBUGPRINT[5]: main.cpp:422: mp=" << mp / g_system.system_int_size() << std::endl;
    auto ds = tree.m_count_ds_inter;
    std::cerr << "DEBUGPRINT[6]: main.cpp:424: ds=" << ds / g_system.system_int_size() << std::endl;
}

auto main(const int argc, const char *const argv[]) -> int {
    if (argc != 2) {
        Logging::err("Please supply a single file argument!");
        return -1;
    }
    // initialize the g_system variable
    g_system.init_system(argv[1]);
    // g_system.init_energy();

    // ============================================================================================
    // task 1
    // ============================================================================================
    // plot_rho_step_1();
    // TODO: (aver)
    // - We still need to do a comparison between different softeining values and discuss their
    //      significance
    // - Also explain dependence of force calculation on direct force calculation
    plot_forces_step_2();

    g_system.calc_real_relaxation();
    // ============================================================================================
    // task 2
    // ============================================================================================

    // plot_do_steps();

    // Timer tree_timer = Timer();
    // tree_timer.start();
    tree_code();
    // tree_timer.stop();
    // tree_timer.print_duration();
    // plot_gif_steps();
    // g_system.animate_particles();

    // auto etas = {0.1, 0.01};
    // auto etas = {0.1};
    // auto etas = {0.01};

    // auto divs = {1.};
    // auto divs = {100.};
    // auto divs = {250,};
    // for (auto eta : etas) {
    //     for (auto div : divs) {
    //         auto direct_en = plot_gif_steps(eta, div);
    //         plot_en(std::get<0>(direct_en), std::get<1>(direct_en), eta, div);
    //         plot_en(std::get<0>(direct_en), std::get<2>(direct_en), eta, div, true);
    //     }
    // }
    // auto tree_en = plot_gif_steps_tree();
    // comp_cost();

    Logging::info("Successfully quit!");
    return 0;
}
