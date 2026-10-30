#pragma once

#include <CGAL/Simplification_2/data_structures/point_quad_tree.h>

namespace CGAL::internal {

	template<class Graph> 
	struct VertexQuadTreeTraits {

		using Element = Graph::Vertex_handle;
		using Kernel = Graph::Kernel;

		static const Point<Kernel>& get_point(Element elt) {
			return elt->point();
		}
	};	

	template<class Graph>
	using VertexQuadTree = CGAL::internal::PointQuadTree<VertexQuadTreeTraits<Graph>>;
}