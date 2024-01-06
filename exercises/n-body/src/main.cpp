#include "data.hpp"
#include "logging.hpp"
#include "particle.hpp"
#include "types.hpp"

#include <mgl2/mgl.h>

#include <future>
#include <ranges>

static PartVec particles;

auto plot_part(mglGraph *gr) {
    auto future_x = std::async(std::launch::async, [&]() {
        auto transformed_x = particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.x();
                             });
        std::vector<double> transformed_vector_x(transformed_x.begin(), transformed_x.end());
        return transformed_vector_x;
    });
    auto future_y = std::async(std::launch::async, [&]() {
        auto transformed_y = particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.y();
                             });
        std::vector<double> transformed_vector_y(transformed_y.begin(), transformed_y.end());
        return transformed_vector_y;
    });
    auto future_z = std::async(std::launch::async, [&]() {
        auto transformed_z = particles | std::views::transform([&](const Particle3D &particle) {
                                 return particle.position.z();
                             });
        std::vector<double> transformed_vector_z(transformed_z.begin(), transformed_z.end());
        return transformed_vector_z;
    });

    mglData x = future_x.get();
    mglData y = future_y.get();
    mglData z = future_z.get();

    auto furthest_particle = Particles::get_max_distance(particles);

    // set plot parameters
    gr->SetSize(1920, 1080);
    gr->Rotate(50, 60);
    gr->SetRanges(furthest_particle.position.x(),
                  -furthest_particle.position.x(),
                  furthest_particle.position.y(),
                  -furthest_particle.position.y(),
                  furthest_particle.position.z(),
                  -furthest_particle.position.z());
    gr->Box();
    gr->Axis("xyz");
    gr->Alpha(true);
    // Draw the sphere at the origin with the specified radius
    gr->Sphere(mglPoint(0, 0, 0), furthest_particle.distance, "b");
    gr->Sphere(mglPoint(0, 0, 0), furthest_particle.distance / 2, "g");
    gr->Alpha(false);
    gr->Dots(x, y, z, "r");
    gr->WriteFrame("plot.png");
    return 0;
}

auto main(int argc, char *argv[]) -> int {
    if (argc != 2) {
        Logging::err("File argument missing");
        return -1;
    }

    auto particles_opt = Data::read_data(argv[1]);
    if (not particles_opt.has_value()) {
        Logging::err("Error while reading file");
        return -1;
    };

    // auto particles = particles_opt.value();
    particles = particles_opt.value();
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
