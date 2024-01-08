#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "system.hpp"

#include <cmath>
#include <format>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

auto plot_part(mglGraph *gr) {
    // set plot parameters
    gr->SetSize(1920, 1080);

    mglData x;
    mglData y;
    mglData z;

    auto trans_vec = g_system.transform_vectors();

    x = std::get<0>(trans_vec);
    y = std::get<1>(trans_vec);
    z = std::get<2>(trans_vec);

    auto furthest_particle = g_system.get_max_distance();

    // auto x_bound = std::abs(furthest_particle.position.x());
    // auto y_bound = std::abs(furthest_particle.position.y());
    // auto z_bound = std::abs(furthest_particle.position.z());

    Histogram hist(1'000'000, furthest_particle.distance, g_system);

    g_system.update_half_mass(hist.m_shells);
    g_system.update_scale_length();

    // set plot parameters
    gr->SetSize(1920, 1080);
    Logging::info("Total mass of system: {}", g_system.m_total_mass);
    Logging::info("Half mass of system: {}", g_system.m_half_mass);
    Logging::info("Scaling length of system: {}", g_system.m_scale_length);

    std::vector<double> index;
    std::vector<double> hernquist_dens;
    std::vector<double> numeric_dens;
    constexpr auto no_steps = 50;
    const auto k_uniform_mass = g_system.m_particles[1].mass;
    g_system.m_min_rad = 0.005;

    auto drLinToLog = [](const double i) {
        return g_system.m_min_rad * std::pow(g_system.m_max_rad / g_system.m_min_rad, i / no_steps);
    };
    auto plotFit = [](double x) { return x + std::numeric_limits<double>::epsilon(); };

    // mglData y2(no_steps);

    // precalculated (by the compiler!) constant for shell volume
    constexpr auto K_SHELL_VOL_PREF = 4. / 3. * std::numbers::pi;

    for (int r = 0; r <= no_steps; r++) {
        auto lower_rad = drLinToLog(r);
        auto upper_rad = drLinToLog(r + 1);

        // auto lower_rad = static_cast<double>(r);
        // auto upper_rad = static_cast<double>(r + 1);

        std::cerr << "DEBUGPRINT[7]: main.cpp:112: lower_rad=" << lower_rad << std::endl;
        std::cerr << "DEBUGPRINT[6]: main.cpp:113: upper_rad=" << upper_rad << std::endl;

        // y2.a[r] = plotFit(g_system.density_hernquist((lower_rad / (upper_rad + 1)) * .5));

        const auto k_shell_mass =
            g_system.get_constrained_mass(upper_rad) - g_system.get_constrained_mass(lower_rad);
        // std::cerr << "DEBUGPRINT[5]: main.cpp:121: k_shell_mass=" << k_shell_mass << std::endl;

        const auto k_shell_volume =
            K_SHELL_VOL_PREF * (std::pow(upper_rad, 3) - std::pow(lower_rad, 3));
        // std::cerr << "DEBUGPRINT[4]: main.cpp:124: k_shell_volume=" << k_shell_volume <<
        // std::endl;

        const auto k_shell_rho = k_shell_mass / k_shell_volume;
        // std::cerr << "DEBUGPRINT[3]: main.cpp:127: k_shell_rho=" << k_shell_rho << std::endl;

        const auto k_no_parts_in_shell = k_shell_mass / k_uniform_mass;
        // std::cerr << "DEBUGPRINT[2]: main.cpp:129: k_no_parts_in_shell=" << k_no_parts_in_shell
        //           << std::endl;

        const auto k_rho_error = std::sqrt(k_no_parts_in_shell) * k_uniform_mass / k_shell_volume;
        // std::cerr << "DEBUGPRINT[1]: main.cpp:131: k_rho_error=" << k_rho_error << std::endl;

        // y2.a[r] = g_system.density_hernquist((lower_rad + upper_rad) / 2);
        // Logging::dbg(y2.a[r]);

        auto val = plotFit(g_system.density_hernquist((lower_rad + upper_rad) * .5));
        std::cerr << "DEBUGPRINT[2]: main.cpp:139: val=" << val << std::endl;

        index.emplace_back(r);
        hernquist_dens.emplace_back(plotFit(val));
        numeric_dens.emplace_back(plotFit(k_shell_rho));

        Logging::info("\n");
    }

    mglData x2 = index;
    mglData y2 = hernquist_dens;
    mglData y3 = numeric_dens;
    Logging::info("{}", x2.Maximal());

    // gr->SetRanges(y2.Minimal(), y2.Maximal());
    // gr->SetRange('y', y2.Minimal(), y2.Minimal());
    // gr->SetRange('x', x2);

    auto outMin = std::min(y2.Minimal(), y3.Minimal());
    auto outMax = std::max(y2.Maximal(), y3.Maximal());
    gr->SetRange('x', x2);
    gr->SetRange('y', outMin, outMax);

    gr->Axis();
    gr->Plot(x2, y2, "b");
    gr->AddLegend("Hernquist Density Profile", "b");
    gr->Plot(x2, y3, "r. ");
    gr->AddLegend("Numeric Density Profile", "r");
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

#if 0
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
