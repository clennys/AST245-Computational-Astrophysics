#include "treecode.hpp"
#include "logging.hpp"

#include <cassert>

TreeCode::TreeCode(BoundingCube root_box, PartVec particles, double tolerance_angle)
    : m_bounding_cube_root(root_box), m_particles(particles), m_tolerance_angle(tolerance_angle) {
    m_octree = new Node(nullptr, m_particles, m_bounding_cube_root, 0);
}

TreeCode::~TreeCode() { delete m_octree; }

auto TreeCode::recursive_build_tree(Node *node) -> void {
    m_height = std::max(node->populate_children(), m_height);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = node->m_children(i, j, k);
                if (child && child->m_particles.size() > 1) {
                    recursive_build_tree(child);
                }
            }
        }
    }
}

auto TreeCode::build() -> void {
    Logging::info("Start Building Octree...");
    recursive_build_tree(m_octree);
    Logging::info("Done with height {}", m_height);
}

auto TreeCode::recursive_tree_walk(Node *node, const Particle3D &part) -> Eigen::Vector3d {
    auto opening_angle = node->calc_opening_angle(part);

    // TODO: (aver) Properly check for a leaf node, in this case we assume a leaf node consists of a
    // singular particle in the node
    if (node->m_particles.size() < 2) {
        return part.calc_direct_force_with_part(node->m_particles[0]);
    }

    if (opening_angle < m_tolerance_angle) {
        return node->multipole_expansion(part);
    } else {
        // Logging::warn("Opening angle ignored with {} particles", node->m_particles.size());
        // WARN: (aver) We shouldn't be able to fall into this case
        assert(node->m_particles.size() > 1);
    }

    auto force = Eigen::Vector3d().setZero();

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Node *child = node->m_children(i, j, k);
                if (child) {
                    force += recursive_tree_walk(child, part);
                }
            }
        }
    }

    return force;
}

auto TreeCode::tree_walk() -> void {
    // TODO: (aver) add parallelisation, consider, that Node is a complex type, with a parent
    // pointer, simply copying this->m_octree to a new pointer results in a segfault
    for (auto &part : m_particles) {
        part.m_tree_force = recursive_tree_walk(m_octree, part);
    }
}

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

auto TreeCode::plot(mglGraph &gr) -> void { plot_recursive(gr, m_octree); };
