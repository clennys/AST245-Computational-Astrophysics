#include "particle.hpp"
#include "logging.hpp"

#include <iostream>

auto Particle3D::calc_origin_dist() -> double { return this->position.norm(); }

auto Particle3D::update_origin_dist() -> void { this->distance = calc_origin_dist(); }

auto Particle3D::print_summary() const -> void {
    Logging::info("Printing summary:");
    std::cout << "Position =\n" << this->position << std::endl;
    std::cout << "Velocity =\n" << this->velocity << std::endl;
    std::cout << "Potential = " << this->potential << std::endl;
    std::cout << "Distance = " << this->distance << std::endl;
    Logging::info("Summary DONE");
}
