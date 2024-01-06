#include "data.hpp"
#include "histogram.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "types.hpp"

#include <cstdlib>
#include <format>
#include <mgl2/mgl.h>

#include <future>
#include <ranges>

/// Particles read in and used in the tasks.
/// Global because the MathGL functions are not allowed to take parameters
///
/// TODO: (aver) consider hiding this/reducing its scope
static PartVec g_particles;

auto plot_part(mglGraph *gr) {
    // TODO: (aver) consider moving the transformations into a separate function
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

    Histogram hist(30, furthest_particle.distance, g_particles);

    // set plot parameters
    gr->SetSize(1920, 1080);
    // gr->Rotate(50, 60);
    gr->SetRanges(-x_bound, x_bound, -y_bound, y_bound, -z_bound, z_bound);
    gr->Box();
    gr->Axis("xyz");
    gr->Alpha(true);

    // Draw the sphere at the origin with the specified radius
    std::vector<double> v_bin_value, v_bin_idx;
    {
        // gr->Sphere(mglPoint(0, 0, 0), furthest_particle.distance, "b");
        // gr->Sphere(mglPoint(0, 0, 0), furthest_particle.distance / 2, "g");
        const mglPoint origin(0, 0, 0);

        for (auto &shell : hist.m_shells) {
            // Logging::info(shell.m_particles.size());
            gr->Sphere(origin, shell.m_upper, "b");

            Logging::dbg(std::format("val: {}, idx: {}", shell.m_particles.size(), shell.m_index));

            // populate values for plot
            v_bin_value.emplace_back(static_cast<int>(shell.m_particles.size()));
            v_bin_idx.emplace_back(static_cast<int>(shell.m_index));
        }
    }
    // convert values to mgl types
    mglData bin_value = v_bin_value;
    mglData bin_idx = v_bin_idx;

    gr->Alpha(false);
    gr->Dots(x, y, z, "r");
    gr->WriteFrame("plot.png");

    gr->ClearFrame();
    gr->ResetFrames();

    gr->SetRanges(bin_idx.Minimal(), bin_idx.Maximal(), bin_value.Minimal(), bin_value.Maximal());
    gr->Axis("xy");
    gr->Bars(bin_idx, bin_value);
    gr->WriteFrame("histogram.png");
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

    // auto particles = particles_opt.value();
    g_particles = particles_opt.value();
    mglGraph gr;
    plot_part(&gr);

    // mglGraph gr(plot_part, "MathGL examples");

    // auto gr_res = gr.Run();
    // if (gr_res != 0) {
    //   Logging::err("MathGL failed");
    //   return gr_res;
    // }
    //
    // Logging::info("Successfully quit!");
    // Logging::info(Particles::get_max_distance(particles));
    return 0;
}
