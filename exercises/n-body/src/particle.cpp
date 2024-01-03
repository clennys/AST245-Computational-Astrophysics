#include "particle.hpp"
#include "logging.hpp"

#include <algorithm>

auto Particle3D::calc_orign_distance() -> void {
    this->distance = this->position.norm();
}

auto Particles::get_max_distance(const std::vector<Particle3D> &particles) -> Particle3D {
    auto result =
        std::max_element(particles.begin(), particles.end(), [](Particle3D a, Particle3D b) {
            return a.distance < b.distance;
        });
    Logging::info(result->position);
    Logging::info(result->distance);
    return *result;
}
