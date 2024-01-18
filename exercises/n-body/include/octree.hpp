#ifndef OCTREE_H_
#define OCTREE_H_

#include "node.hpp"

class Octree {
	public:
		Node* m_root_node;
		BoundingCube m_bounding_cube_root;
		PartVec m_particles;
		// int m_height = -1;


		Octree(BoundingCube root_box, PartVec particles);
		auto build() -> void;
		auto recursive_build_tree(Node *root) -> void;

};


#endif // ! OCTREE_H_
