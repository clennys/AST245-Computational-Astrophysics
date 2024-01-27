#ifndef NODE_H_
#define NODE_H_

#include "particle.hpp"

#include <unsupported/Eigen/CXX11/Tensor>

class Node;

using BoundingCube = Eigen::Tensor<Eigen::Vector3d, 3>;
using SubBoundingCubes = Eigen::Tensor<BoundingCube, 3>;
using SubCubesNodes = Eigen::Tensor<Node *, 3>;

class Node {
  public:
    // parent pointer
    Node *m_parent = nullptr;
    // array of children (8)
    SubCubesNodes m_children;

    // particles in node
    PartVec m_particles = {};

    // Bounding cub index (x,y,z)
    //    (0,1,1) .+------+ (1,1,1)
    //          .' |    .'|
    // (0,1,0) +------+' (1,1,0)
    //         |   |  |   | y
    //     (0,0,1).+--+---+ (1,0,1)   | z
    //         |.'    | .'|/
    // (0,0,0) +------+' (1,0,0) +--x
    //
    BoundingCube m_bounding_cube;
    int m_depth = -1;
    Eigen::Matrix3d m_quadrupole;
    double m_monopole = 0;
    Eigen::Vector3d m_center_of_mass;

    explicit Node(Node *par, PartVec part_vec, BoundingCube cube, int depth);
    ~Node();

    Node() = delete;
    Node(const Node &) = delete;
    Node &operator=(const Node &) = delete;
    Node(Node &&) = delete;
    Node &operator=(Node &&) = delete;

    auto in_bounding_box(Particle3D part) -> bool;
    auto multipole_expansion(const Particle3D &part) -> Eigen::Vector3d;
    auto calc_expansion_factors() -> void;
    auto calc_opening_angle(const Particle3D &part) const -> double;
    auto populate_children() -> int;
    auto octa_split_bounding_box() -> SubBoundingCubes;
    auto create_sub_bounding_cube(Eigen::Vector3d origin, double cube_side_length) -> BoundingCube;
    auto in_bounding_box(const BoundingCube cube, const Particle3D &part) -> bool;
    auto print_cube() -> void;
};

#endif // ! NODE_H_