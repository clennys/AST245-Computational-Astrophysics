#include "data.hpp"
#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "types.hpp"

#include <cstdlib>
#include <format>
#include <mgl2/mgl.h>
#include <mgl2/qt.h>

#include <future>
#include <numbers>
#include <ranges>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
/// TODO: (aver) consider hiding this/reducing its scope
static PartVec g_particles;

auto density_hernquist(double rad, double total_mass, double scale_length) -> double {
    return (total_mass / (2 * std::numbers::pi)) * (scale_length / rad) *
           (1 / std::pow(rad + scale_length, 3));
}

auto plot_part(mglGraph *gr) {
    // TODO: (aver) consider moving the transformations into a separate function to reduce
    // recalculation on each viewport change
    auto future_x = std::async(std::launch::async, [&]() {
        auto transformed_x = g_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.x();
                             });
        std::vector<double> transformed_vector_x(transformed_x.begin(), transformed_x.end());
        return transformed_vector_x;
    });
    auto future_y = std::async(std::launch::async, [&]() {
        auto transformed_y = g_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.y();
                             });
        std::vector<double> transformed_vector_y(transformed_y.begin(), transformed_y.end());
        return transformed_vector_y;
    });
    auto future_z = std::async(std::launch::async, [&]() {
        auto transformed_z = g_particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.z();
                             });
        std::vector<double> transformed_vector_z(transformed_z.begin(), transformed_z.end());
        return transformed_vector_z;
    });

    mglData x = future_x.get();
    mglData y = future_y.get();
    mglData z = future_z.get();

    auto furthest_particle = Particles::get_max_distance(g_particles);

    auto x_bound = std::abs(furthest_particle.position.x());
    auto y_bound = std::abs(furthest_particle.position.y());
    auto z_bound = std::abs(furthest_particle.position.z());

    Histogram hist(100'000, furthest_particle.distance, g_particles);
    auto half_mass_rad = hist.calc_half_mass();
    auto scale_length_a = half_mass_rad / (1 + (std::sqrt(2)));

    // auto hernquist_dens_profile = density_hernquist()

    Logging::info(std::format("Total mass of system: {}", Particles::g_total_mass));
    Logging::info(std::format("Half mass of system: {}", half_mass_rad));
    Logging::info(std::format("Scaling length of system: {}", scale_length_a));

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

    auto particles_opt = Data::read_data(argv[1]);
    if (not particles_opt.has_value()) {
        Logging::err(std::format("Error while reading file: {}", argv[1]));
        return -1;
    };

    g_particles = particles_opt.value();

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
