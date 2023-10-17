#include "ode_solver.hpp"
#include "particle.hpp"

#include "matplotlibcpp.h"

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <iostream>
#include <ranges>
#include <string>
#include <vector>

// Examples for Eigen
// https://bitbucket.org/chris_richardson/eigen-demo/src/master/explicit.cpp

namespace plt = matplotlibcpp;

// TODO: (aver) add kepler function and fix usage of acceleration
auto main() -> int {
    const auto eccentricities = {0., .8};
    const auto timesteps = {0.001, 0.01};

    ODESolver ode_system(ODESolver::ODEScheme::ExplicitEuler,
                         ODESolver::ODEDerivatives::KeplerianOrbits, 0.01);

    for (auto eccentricity : eccentricities) {
        Particle init_particle(
            Eigen::Vector2d{1, 0},
            Eigen::Vector2d({0, std::sqrt(1 + eccentricity)}), 1.);

        auto particles = ode_system.solve_system(init_particle, 10);

        // for (auto particle : particles) {
        //     std::cout << "position:" << particle.position << "\n";
        //     std::cout << "velocity:" << particle.velocity << '\n';
        // }

        // we need to transform the vectors of particles into positions for
        // plotting
        auto transformed_x =
            particles | std::views::transform([](const Particle &particle) {
                return particle.position.x();
            });
        std::vector<double> transformed_vector_x(transformed_x.begin(),
                                                 transformed_x.end());

        auto transformed_y =
            particles | std::views::transform([](const Particle &particle) {
                return particle.position.y();
            });
        std::vector<double> transformed_vector_y(transformed_y.begin(),
                                                 transformed_y.end());

        plt::plot(transformed_vector_x, transformed_vector_y,
                  std::map<std::string, std::string>{{"marker", "."}});

        plt::scatter(std::vector{0.}, std::vector{0.});
        plt::show();
    }

    // plt::scatter(std::vector{0.}, std::vector{0.});
    // plt::plot( std::vector{ 0., 3., }, std::vector{0., 3.});
    // plt::show();
}
