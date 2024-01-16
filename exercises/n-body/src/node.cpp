#include "node.hpp"
#include <iostream>
#include <memory>

Node::Node() {}

Node::Node(Node *par, PartVec part_vec, BoundingCube box)
    : m_parent(par), m_particles(part_vec), m_bounding_cube(box) {
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
auto Node::in_bounding_box(const BoundingCube cube, const Particle3D part) -> bool {
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

auto Node::populate_children() -> void {
    SubBoundingCubes new_bounding_cubes = octa_split_bounding_box();

    for (const auto &part : m_particles) {
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    if (in_bounding_box(new_bounding_cubes(i, j, k), part)) {
                        if (m_children(i, j, k) != nullptr) {
                            PartVec child_part = {part};
                            m_children(i, j, k) =
                                new Node(this, child_part, new_bounding_cubes(i, j, k));
                        } else {
                            m_children(i, j, k)->m_particles.push_back(part);
                        }
                    }
                }
            }
        }
    }
}
