#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "system.hpp"

#include <cmath>
#include <cstdlib>
#include <format>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
static System g_system;

auto plot_part(mglGraph *gr) {
    auto trans_vec = g_system.transform_vectors();

    mglData x = std::get<0>(trans_vec);
    mglData y = std::get<1>(trans_vec);
    mglData z = std::get<2>(trans_vec);

    auto furthest_particle = g_system.get_max_distance();

    auto x_bound = std::abs(furthest_particle.position.x());
    auto y_bound = std::abs(furthest_particle.position.y());
    auto z_bound = std::abs(furthest_particle.position.z());

    Histogram hist(1000, furthest_particle.distance, g_system);

    // auto half_mass_rad = g_system.calc_half_mass(hist.m_shells);
    // auto scale_length_a = half_mass_rad / (1 + (std::sqrt(2)));
    g_system.update_half_mass(hist.m_shells);
    g_system.update_scale_length();

    // auto hernquist_dens_profile = density_hernquist()

    Logging::info(std::format("Total mass of system: {}", g_system.m_total_mass));
    Logging::info(std::format("Half mass of system: {}", g_system.m_half_mass));
    Logging::info(std::format("Scaling length of system: {}", g_system.m_scale_length));

    // set plot parameters
    gr->SetSize(1920, 1080);

    // set sphere 3d plot params
    // TODO: (aver) fix plot rotation in subsequent plots
    // or just use 2 subplots
    // gr->Rotate(50, 60);
    gr->SetRanges(-x_bound, x_bound, -y_bound, y_bound, -z_bound, z_bound);
    gr->Box();
    gr->Axis("xyz");
    gr->Alpha(true);

    Logging::info("Plotting Spheres and preparing histogram...");

    // Draw the sphere at the origin with the specified radius
    std::vector<double> v_bin_value, v_bin_idx;
    { // reduced scope, so that origin is invalidated afterwards
        const mglPoint origin(0, 0, 0);

        for (auto &shell : hist.m_shells) {
            // plot bin sphere
            // TODO: (aver) adjust color and transpency, otherwise not much sense
            gr->Sphere(origin, shell.m_upper, "b");

            // Logging::dbg(std::format("val: {}, idx: {}", shell.m_particles.size(),
            // shell.m_index));

            // populate values for plot
            v_bin_value.emplace_back(static_cast<int>(shell.m_particles.size()));
            v_bin_idx.emplace_back(static_cast<int>(shell.m_upper));
        }
    }
    gr->Alpha(false);
    gr->Dots(x, y, z, "r");
    gr->WriteFrame("plot.png");
    Logging::info("Spheres plotted.");

    // reset frames and set options for histogram plot
    gr->ClearFrame();
    gr->ResetFrames();

    // convert values to mgl types
    x = v_bin_idx;
    y = v_bin_value;

    gr->SetRanges(x.Minimal(), x.Maximal(), y.Minimal(), y.Maximal());
    gr->Axis("xy");
    gr->Bars(x, y);
    gr->WriteFrame("histogram.png");
    Logging::info("Histogram plotted.");
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
