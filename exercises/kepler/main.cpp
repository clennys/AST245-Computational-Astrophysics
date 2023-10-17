#include "ode_solver.hpp"
#include "particle.hpp"

#include "matplotlibcpp.h"

#include <Eigen/Dense>
#include <cmath>
#include <functional>
#include <iostream>
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
            Eigen::Vector2d{0, 1},
            Eigen::Vector2d({0, std::sqrt(1 + eccentricity)}), 1.);

        auto particles = ode_system.solve_system(init_particle, 10);
        std::cerr << "DEBUGPRINT[1]: kepler.cpp:32: particles="
                  << particles.size() << std::endl;

        for (auto particle : particles) {
            std::cout << particle.position << '\n';
            std::cout << particle.velocity << '\n';
        }
    }
    plt::scatter(std::vector{0.}, std::vector{0.});
    plt::show();
}
