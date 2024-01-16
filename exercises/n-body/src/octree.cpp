#include "octree.hpp"

Octree::Octree(BoundingCube root_box, PartVec particles)
    : m_bounding_cube_root(root_box), m_particles(particles) {
    m_root_node = Node(nullptr, particles, root_box);
}

auto Octree::build() -> void {}


