#ifndef TREECODE_H_
#define TREECODE_H_

#include "node.hpp"
#include <mgl2/mgl.h>

class TreeCode {
  public:
    Node *m_octree;
    BoundingCube m_bounding_cube_root;
    PartVec m_particles;
    int m_height = -1;
    double m_tolerance_angle;

    explicit TreeCode(BoundingCube root_cube, PartVec particles, double tolerance_angle);
    ~TreeCode();
    auto build() -> void;
    auto recursive_build_tree(Node *root) -> void;
    auto tree_walk() -> void;
    auto recursive_tree_walk(Node *root, const Particle3D &part) -> Eigen::Vector3d;
    auto plot(mglGraph &gr) -> void;
    auto plot_recursive(mglGraph &gr, const Node *node) -> void;
    auto plot_cube(mglGraph &gr, const Node *node) -> void;
};

#endif // ! TREECODE_H_
