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
		int m_count_mp_inter = 0;
		int m_count_ds_inter = 0;
		double m_pot_energy = 0.;
		double m_kin_energy = 0.;
		double m_tot_energy = 0.;

    explicit TreeCode(BoundingCube root_cube, PartVec particles, double tolerance_angle);
    ~TreeCode();
    auto build() -> void;
    auto recursive_build_tree(Node *root) -> void;
    auto tree_walk() -> void;
    auto tree_step(const double dt) -> void;
    auto recursive_tree_walk(Node *root, const Particle3D &part) -> Eigen::Vector3d;
    auto plot(mglGraph &gr) -> void;
    auto plot_recursive(mglGraph &gr, const Node *node) -> void;
    auto plot_cube(mglGraph &gr, const Node *node) -> void;
    auto reset_tree() -> void;
    auto computational_cost(double radius, double factor) -> double;
		auto total_energy(double k_dim_mass) -> void;
};

#endif // ! TREECODE_H_
