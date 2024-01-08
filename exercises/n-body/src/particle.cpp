#include "particle.hpp"

#include <iostream>

auto Particle3D::calc_orign_distance() -> void { this->distance = this->position.norm(); }

auto Particle3D::print_summary() const -> void {
    std::cout << "Position = " << this->position << std::endl;
    std::cout << "Potential = " << this->potential << std::endl;
    std::cout << "Velocity = " << this->velocity << std::endl;
    std::cout << "Distance = " << this->distance << std::endl;
}
