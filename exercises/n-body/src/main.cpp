#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "system.hpp"

#include <cmath>
#include <format>
#include <mgl2/mgl.h>
// #include <mgl2/qt.h>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

auto plot_part(mglGraph *gr) {
    // set plot parameters
    gr->SetSize(1920, 1080);

    auto furthest_particle = g_system.get_max_distance();
    Histogram hist(100'000, furthest_particle.distance, g_system);

    // g_system.update_half_mass(hist.m_shells);
    g_system.update_half_mass_radius(hist.m_shells);
    g_system.update_scale_length();

    // set plot parameters
    gr->SetSize(1920, 1080);
    Logging::info("Total mass of system: {}", g_system.m_total_mass);
    Logging::info("Half mass of system: {}", g_system.m_half_mass_rad);
    Logging::info("Scaling length of system: {}", g_system.m_scale_length);

    std::vector<double> index;
    std::vector<double> hernquist_dens;
    std::vector<double> numeric_dens;

    // precalculated (by the compiler!) constant for shell volume
    constexpr auto k_shell_vol_pref = 4. / 3. * std::numbers::pi;
#if 0
    auto i = 0;
    for (const auto &shell : hist.m_shells) {

        const auto k_shell_volume =
            k_shell_vol_pref * (std::pow(shell.m_upper, 3) - std::pow(shell.m_lower_inc, 3));

        const auto k_shell_rho = shell.m_mass / k_shell_volume;

        auto val = g_system.density_hernquist((shell.m_lower_inc + shell.m_upper) * .5);

        Logging::info("\n\tH: {}\n\tN: {}", val, k_shell_rho);

        index.emplace_back(i);
        hernquist_dens.emplace_back(val);
        numeric_dens.emplace_back(k_shell_rho);
        i++;
    }
#else
    constexpr auto no_steps = 50;
    g_system.m_min_rad = 0.005;

    // TODO: (aver) make these static/member methods of gsystem
    constexpr auto dr_lin_to_log = [&](const double i) {
        return g_system.m_min_rad * std::pow(g_system.m_max_rad / g_system.m_min_rad, i / no_steps);
    };
    constexpr auto fit_for_plot = [](const double x) {
        return x + std::numeric_limits<double>::epsilon();
    };

    for (int r = 0; r <= no_steps; r++) {
        auto lower_rad = dr_lin_to_log(r);
        auto upper_rad = dr_lin_to_log(r + 1);
        const auto k_shell_mass = g_system.get_constrained_shell_mass(lower_rad, upper_rad);
        const auto k_shell_volume =
            k_shell_vol_pref * (std::pow(upper_rad, 3) - std::pow(lower_rad, 3));
        const auto k_no_parts_in_shell = k_shell_mass / g_system.km_mass;
        const auto k_rho_error = std::sqrt(k_no_parts_in_shell) * g_system.km_mass / k_shell_volume;

        const auto k_shell_rho = fit_for_plot(k_shell_mass / k_shell_volume);
        const auto k_hern_rho =
            fit_for_plot(g_system.density_hernquist((lower_rad + upper_rad) / 2));

        Logging::info("\n\tH: {}\n\tN: {}", k_hern_rho, k_shell_rho);
        index.emplace_back(r);
        hernquist_dens.emplace_back(k_hern_rho);
        numeric_dens.emplace_back(k_shell_rho);
    }
#endif

    mglData x = index;
    mglData y_hern = hernquist_dens;
    mglData y_num = numeric_dens;

    auto y_min = std::min(y_hern.Minimal(), y_num.Minimal());
    auto y_max = std::max(y_hern.Maximal(), y_num.Maximal());

    gr->SetRange('x', x);
    gr->SetRange('y', y_min, y_max);

    gr->Axis();

    gr->Plot(x, y_hern, "b");
    gr->AddLegend("Hernquist Density Profile", "b");

    gr->Plot(x, y_num, "r");
    gr->AddLegend("Numeric Density Profile", "r");

    gr->Legend();
    gr->WriteFrame("hernquist.jpg");
    gr->WriteFrame("hernquist.png");
    Logging::info("Hernquist plotted.");

    return 0;
}

auto main(int argc, char *argv[]) -> int {
    if (argc != 2) {
        Logging::err("File argument missing");
        return -1;
    }

    g_system = System(argv[1]);

#if 1
    mglGraph gr;
    plot_part(&gr);
#else
    mglQT gr(plot_part, "MathGL examples");
    auto gr_res = gr.Run();
    if (gr_res != 0) {
        Logging::err("MathGL failed");
        return gr_res;
    }
#endif

    Logging::info("Successfully quit!");
    return 0;
}
