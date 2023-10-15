#include "particle.hpp"

Particle::Particle(Eigen::Vector2d pos, Eigen::Vector2d vel, double mass)
    : position(pos), velocity(vel), k_mass(mass) {}

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
