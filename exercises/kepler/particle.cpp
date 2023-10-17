#include "particle.hpp"
#include <algorithm>

Particle::Particle(Eigen::Vector2d pos, Eigen::Vector2d vel, double mass)
    : position(pos), velocity(vel), k_mass(mass) {}

Particle &Particle::operator=(const Particle &other) {
    if (this == &other) {
        return *this;
    }
    this->position = other.position;
    this->velocity = other.velocity;
    this->time = other.time;
    this->energy = other.energy;
    // this->k_mass = other.k_mass;
    this->eccentricity = other.eccentricity;
    this->angular_momentum = other.angular_momentum;
    return *this;
}
Particle::~Particle() {}

auto Particle::compute_energy() -> void {
    const auto subtractee = std::pow(velocity.norm(), 2) / 2;
    const auto subtrahend = 1 / position.norm();
    this->energy = subtractee - subtrahend;
}

auto Particle::compute_angular_momentum() -> void {
    const auto left = position.x() * velocity.y();
    const auto right = position.y() * velocity.x();
    this->angular_momentum = left - right;
}

auto Particle::compute_eccentricity() -> void {
    const auto left =
        position * (std::pow(velocity.norm(), 2) - (1 / position.norm()));
    const auto right = velocity * velocity.dot(position);
    this->eccentricity = (left - right).norm();
}
