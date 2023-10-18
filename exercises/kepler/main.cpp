#include "ode_solver.hpp"
#include "particle.hpp"

#include "matplotlibcpp.h"

#include <Eigen/Dense>
#include <cmath>
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

        auto transformed_positions = ODESolver::transform_vec2d(
            particles, ODESolver::TransElemem::Position);


        plt::plot(transformed_positions.first, transformed_positions.second,

                  std::map<std::string, std::string>{{"marker", "."}});

        plt::scatter(std::vector{0.}, std::vector{0.});
        plt::show();
    }
}
