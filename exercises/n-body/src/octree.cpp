#include "octree.hpp"
#include "logging.hpp"
#include <iostream>

Octree::Octree(BoundingCube root_box, PartVec particles)
    : m_bounding_cube_root(root_box), m_particles(particles) {
    m_root_node = new Node(nullptr, m_particles, m_bounding_cube_root, 0);
}

auto Octree::recursive_build_tree(Node *root) -> void {
    root->populate_children();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = root->m_children(i, j, k);
                if (child != nullptr && child->m_particles.size() > 1) {
                    recursive_build_tree(child);
                }
            }
        }
    }
}

auto Octree::build() -> void { recursive_build_tree(m_root_node); }
