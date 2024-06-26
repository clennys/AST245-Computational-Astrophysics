#include "treecode.hpp"
#include "logging.hpp"
#include <numbers>

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
    const auto opening_angle = node->calc_opening_angle(part);

    // TODO: (aver) Properly check for a leaf node, in this case we assume a leaf node consists of a
    // singular particle in the node
    if (node->m_particles.size() < 2) {
        m_count_ds_inter++;
        return part.calc_direct_force_with_part(node->m_particles[0]);
    }

    if (opening_angle < m_tolerance_angle) {
        m_count_mp_inter++;
        // WARN: (dhub) Possible complications if tree-walk is implemented in parallell
        m_force_error += node->force_error(part);
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
    Logging::info("Tree Walk completed");
}

auto TreeCode::tree_step(const double dt) -> void {
    for (auto &part : m_particles) {
        // First leap, get mid velocity
        const auto velocity_mid = part.m_velocity + (part.m_tree_force / part.m_mass) * dt * .5;
        // Update position for force calculation
        part.m_position += velocity_mid * dt;
        // update new force
        part.m_tree_force = recursive_tree_walk(m_octree, part);
        // set new velocity, completing leap-frog
        part.m_velocity = velocity_mid + (part.m_tree_force / part.m_mass) * dt * .5;
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

auto TreeCode::reset_tree() -> void {
    delete this->m_octree;
    this->m_octree = new Node(nullptr, m_particles, m_bounding_cube_root, 0);
}
// NOTE: (dhub) Define the computational cost as the number of nodes
auto TreeCode::computational_cost(double radius, double factor) -> double {
    auto interpart_sep = (((4 * std::numbers::pi) / 3) * radius * radius * radius);
    interpart_sep /= (static_cast<int>(m_particles.size()) * factor);
    interpart_sep = std::pow(interpart_sep, 1. / 3.);
    std::cerr << "DEBUGPRINT[1]: treecode.cpp:158: interpart_sep=" << interpart_sep << std::endl;

    return 4 * std::numbers::pi / std::pow(m_tolerance_angle, 3) * std::log(radius / interpart_sep);
}

auto TreeCode::total_energy(double k_dim_mass) -> void {
    double pot_energy = 0.;
    double kin_energy = 0.;
    for (uint i = 0; i < m_particles.size(); i++) {
        kin_energy += 0.5 * k_dim_mass * m_particles[i].m_velocity.squaredNorm();
        for (uint j = i + 1; j < m_particles.size(); j++) {
            pot_energy += k_dim_mass * k_dim_mass /
                          (m_particles[i].m_position - m_particles[j].m_position).norm();
        }
    }
    m_pot_energy = pot_energy;
    m_kin_energy = kin_energy;
    m_tot_energy = pot_energy - kin_energy;
}
