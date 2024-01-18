#ifndef OCTREE_H_
#define OCTREE_H_

#include "node.hpp"
#include <mgl2/mgl.h>

class Octree {
	public:
		Node* m_root_node;
		BoundingCube m_bounding_cube_root;
		PartVec m_particles;
		// int m_height = -1;


		Octree(BoundingCube root_box, PartVec particles);
		auto build() -> void;
		auto recursive_build_tree(Node *root) -> void;
		auto tree_walk() -> void;
		auto recursive_tree_walk(Node *root) -> void;
		auto plot(mglGraph &gr) -> void;
		auto plot_recursive(mglGraph &gr, const Node *node) -> void;
		auto plot_cube(mglGraph &gr, const Node *node) -> void;

};


#endif // ! OCTREE_H_
