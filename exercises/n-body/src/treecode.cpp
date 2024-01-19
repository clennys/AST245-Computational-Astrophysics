#include "treecode.hpp"
#include "logging.hpp"
#include <iostream>

TreeCode::TreeCode(BoundingCube root_box, PartVec particles, double crit_opening_angle)
    : m_bounding_cube_root(root_box), m_particles(particles), m_crit_opening_angle(crit_opening_angle) {
    m_octree = new Node(nullptr, m_particles, m_bounding_cube_root, 0);
}

auto TreeCode::recursive_build_tree(Node *node) -> void {
    node->populate_children();
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = node->m_children(i, j, k);
                if (child != nullptr && child->m_particles.size() > 1) {
                    recursive_build_tree(child);
                }
            }
        }
    }
}

auto TreeCode::build() -> void { recursive_build_tree(m_octree); }

auto TreeCode::recursive_tree_walk(Node *node) -> void {
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = node->m_children(i, j, k);
                if (child != nullptr) {
                    recursive_tree_walk(child);
                }
            }
        }
    }
}

auto TreeCode::tree_walk() -> void { recursive_tree_walk(m_octree); }

// TODO: (dhub) Use loops to draw lines
auto TreeCode::plot_cube(mglGraph &gr, const Node *node) -> void {
    // Ensure there are exactly 8 corners
    // Draw lines between the vertices to form the cube
    BoundingCube cube = node->m_bounding_cube;
    // draw bottom face
    gr.Line(mglPoint(cube(0, 0, 0).x(), cube(0, 0, 0).y(), cube(0, 0, 0).z()),
            mglPoint(cube(1, 0, 0).x(), cube(1, 0, 0).y(), cube(1, 0, 0).z()));

    gr.Line(mglPoint(cube(1, 0, 0).x(), cube(1, 0, 0).y(), cube(1, 0, 0).z()),
            mglPoint(cube(1, 0, 1).x(), cube(1, 0, 1).y(), cube(1, 0, 1).z()));

    gr.Line(mglPoint(cube(1, 0, 1).x(), cube(1, 0, 1).y(), cube(1, 0, 1).z()),
            mglPoint(cube(0, 0, 1).x(), cube(0, 0, 1).y(), cube(0, 0, 1).z()));

    gr.Line(mglPoint(cube(0, 0, 1).x(), cube(0, 0, 1).y(), cube(0, 0, 1).z()),
            mglPoint(cube(0, 0, 0).x(), cube(0, 0, 0).y(), cube(0, 0, 0).z()));

    // draw top face
    gr.Line(mglPoint(cube(0, 1, 0).x(), cube(0, 1, 0).y(), cube(0, 1, 0).z()),
            mglPoint(cube(1, 1, 0).x(), cube(1, 1, 0).y(), cube(1, 1, 0).z()));

    gr.Line(mglPoint(cube(1, 1, 0).x(), cube(1, 1, 0).y(), cube(1, 1, 0).z()),
            mglPoint(cube(1, 1, 1).x(), cube(1, 1, 1).y(), cube(1, 1, 1).z()));

    gr.Line(mglPoint(cube(1, 1, 1).x(), cube(1, 1, 1).y(), cube(1, 1, 1).z()),
            mglPoint(cube(0, 1, 1).x(), cube(0, 1, 1).y(), cube(0, 1, 1).z()));

    gr.Line(mglPoint(cube(0, 1, 1).x(), cube(0, 1, 1).y(), cube(0, 1, 1).z()),
            mglPoint(cube(0, 1, 0).x(), cube(0, 1, 0).y(), cube(0, 1, 0).z()));

    // draw vertical edges
    gr.Line(mglPoint(cube(0, 0, 0).x(), cube(0, 0, 0).y(), cube(0, 0, 0).z()),
            mglPoint(cube(0, 1, 0).x(), cube(0, 1, 0).y(), cube(0, 1, 0).z()));

    gr.Line(mglPoint(cube(1, 0, 0).x(), cube(1, 0, 0).y(), cube(1, 0, 0).z()),
            mglPoint(cube(1, 1, 0).x(), cube(1, 1, 0).y(), cube(1, 1, 0).z()));

    gr.Line(mglPoint(cube(1, 0, 1).x(), cube(1, 0, 1).y(), cube(1, 0, 1).z()),
            mglPoint(cube(1, 1, 1).x(), cube(1, 1, 1).y(), cube(1, 1, 1).z()));

    gr.Line(mglPoint(cube(0, 0, 1).x(), cube(0, 0, 1).y(), cube(0, 0, 1).z()),
            mglPoint(cube(0, 1, 1).x(), cube(0, 1, 1).y(), cube(0, 1, 1).z()));
}

auto TreeCode::plot_recursive(mglGraph &gr, const Node *node) -> void {
	plot_cube(gr, node);
  for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = node->m_children(i, j, k);
                if (child != nullptr) {
                    plot_recursive(gr, child);
                }
            }
        }
    }
}

auto TreeCode::plot(mglGraph &gr) -> void{
	plot_recursive(gr, m_octree);
};
