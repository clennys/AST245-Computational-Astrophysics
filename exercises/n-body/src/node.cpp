#include "node.hpp"
#include "logging.hpp"
#include "system.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

Node::Node(Node *par, PartVec part_vec, BoundingCube cube, int depth)
    : m_parent(par), m_particles(part_vec), m_bounding_cube(cube), m_depth(depth) {
    m_children = SubCubesNodes(2, 2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                m_children(i, j, k) = nullptr;
            }
        }
    }
}

// TODO: (dhub) Make cleaner.
auto add_to_vect(Eigen::Vector3d vec, double x = 0., double y = 0., double z = 0.)
    -> Eigen::Vector3d {
    Eigen::Vector3d res;
    res << vec.x() + x, vec.y() + y, vec.z() + z;
    return res;
}

// TODO: (dhub) Maybe there's a fancier way, but too tired right now, sry.
auto Node::create_sub_bounding_cube(Eigen::Vector3d origin, double cube_side_length)
    -> BoundingCube {
    BoundingCube cube(2, 2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x_half_side = i == 1 ? cube_side_length : 0.;
                double y_half_side = j == 1 ? cube_side_length : 0.;
                double z_half_side = k == 1 ? cube_side_length : 0.;
                cube(i, j, k) = add_to_vect(origin, x_half_side, y_half_side, z_half_side);
            }
        }
    }
    return cube;
}

auto Node::octa_split_bounding_box() -> SubBoundingCubes {
    SubBoundingCubes sub_boxes(2, 2, 2);

    // NOTE: (dhub) Because we construct cubes all sides have the same length.
    double half_side_length = (m_bounding_cube(1, 0, 0).x() - m_bounding_cube(0, 0, 0).x()) * 0.5;
    // Logging::dbg("Cube sides s halfed: {}", half_side_length);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                double x_half_side = i == 1 ? half_side_length : 0.;
                double y_half_side = j == 1 ? half_side_length : 0.;
                double z_half_side = k == 1 ? half_side_length : 0.;
                Eigen::Vector3d origin =
                    add_to_vect(m_bounding_cube(0, 0, 0), x_half_side, y_half_side, z_half_side);
                sub_boxes(i, j, k) = create_sub_bounding_cube(origin, half_side_length);
            }
        }
    }

    return sub_boxes;
}

// NOTE: (dhub) By construction the cubes are always axis aligned
// Source:
// https://stackoverflow.com/questions/21037241/how-to-determine-a-point-is-inside-or-outside-a-cube
auto Node::in_bounding_box(const BoundingCube cube, const Particle3D &part) -> bool {
    const auto x = part.m_position.x();
    const auto y = part.m_position.y();
    const auto z = part.m_position.z();

    const auto x_min = cube(0, 0, 0).x();
    const auto x_max = cube(1, 0, 0).x();
    const auto y_min = cube(0, 0, 0).y();
    const auto y_max = cube(0, 1, 0).y();
    const auto z_min = cube(0, 0, 0).z();
    const auto z_max = cube(0, 0, 1).z();

    return (x_min <= x && x <= x_max) && (y_min <= y && y <= y_max) && (z_min <= z && z <= z_max);
}

auto Node::populate_children() -> int {
    auto max_depth = 0;
    SubBoundingCubes new_bounding_cubes = octa_split_bounding_box();
    for (const auto &part : m_particles) {
        // NOTE: Particles only fall in single sub bounding cube, there are some on the boundary
        bool particle_allocated = false;
        for (int i = 0; i < 2 && !particle_allocated; i++) {
            for (int j = 0; j < 2 && !particle_allocated; j++) {
                for (int k = 0; k < 2 && !particle_allocated; k++) {
                    if (in_bounding_box(new_bounding_cubes(i, j, k), part)) {
                        if (m_children(i, j, k) == nullptr) {
                            PartVec child_part = {part};
                            m_children(i, j, k) = new Node(
                                this, child_part, new_bounding_cubes(i, j, k), m_depth + 1);
                            max_depth = std::max(max_depth, m_depth + 1);
                        } else {
                            m_children(i, j, k)->m_particles.push_back(part);
                        }
                        particle_allocated = true;
                    }
                }
            }
        }
    }
    return max_depth;
}
auto Node::print_cube() -> void {
    Logging::info("Printing summary");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                Logging::info("Coord Cube at ({},{},{}) is", i, j, k);
                std::cout << m_bounding_cube(i, j, k) << std::endl;
            }
        }
    }
}

auto Node::calc_expansion_factors() -> void {
    static constexpr auto kronecker_delta = [](const int i, const int j) { return i == j ? 1 : 0; };

    m_quadrupole = Eigen::Matrix3d::Zero();
    m_monopole = 0;
    m_center_of_mass = Eigen::Vector3d::Zero();

    for (const auto &part : m_particles) {
        // m_monopole += part.m_mass;
        // m_center_of_mass += part.m_position * part.m_mass;
        m_monopole += System::k_non_dim_mass;
        m_center_of_mass += part.m_position;
    }

    m_center_of_mass /= m_monopole;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            for (const auto &part : m_particles) {
                Eigen::Vector3d d_part_com = m_center_of_mass - part.m_position;
                m_quadrupole(i, j) +=
                    // part.m_mass * (3 * d_part_com(i) * d_part_com(j) -
                    //                kronecker_delta(i, j) * d_part_com.squaredNorm());
                    (3 * d_part_com(i) * d_part_com(j) -
                     kronecker_delta(i, j) * d_part_com.squaredNorm());
            }
        }
    }
}

// TODO: (dhub) Verify formula -> In book they use radius, can we use the cube side length?
// opening angle should be << 1
auto Node::calc_opening_angle(const Particle3D &part) const -> double {
    // NOTE: (dhub) Per construction this is always positive cube(1,0,0).x() >= cube(0,0,0).(x)
    double cube_side = m_bounding_cube(1, 0, 0).x() - m_bounding_cube(0, 0, 0).x();
    double dist_part_com = (part.m_position - m_center_of_mass).norm();
    // auto cube_side / dist_part_com;
    return std::abs(std::atan2(cube_side, dist_part_com));
}

auto Node::multipole_expansion(const Particle3D &part) -> Eigen::Vector3d {
    this->calc_expansion_factors();
    auto dist_p_com = part.m_position - this->m_center_of_mass;
    auto dist_norm = dist_p_com.norm();

    Eigen::Vector3d f_monopole_vec =
        // (m_monopole * dist_p_com * part.m_mass) / (dist_norm * dist_norm * dist_norm);
        -(m_monopole * dist_p_com) / (dist_norm * dist_norm * dist_norm);

    assert(not f_monopole_vec.hasNaN());

    const auto distnorm4inverse = 1 / (dist_norm * dist_norm * dist_norm * dist_norm);
    const auto Qy = (m_quadrupole * dist_p_com);
    const double yQy = dist_p_com.transpose() * Qy;
    const Eigen::Vector3d Q_magn = ((-5. / 2.) * yQy * distnorm4inverse) * dist_p_com;

    const Eigen::Vector3d f_quadrupole_vec = Qy * distnorm4inverse + Q_magn;

    return f_monopole_vec + f_quadrupole_vec;
}

Node::~Node() {
    // Logging::dbg("Destructor of Node called");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            for (int k = 0; k < 2; k++) {
                if (m_children(i, j, k)) {
                    // Logging::dbg("children deleted");
                    delete m_children(i, j, k);
                }
            }
        }
    }
}
