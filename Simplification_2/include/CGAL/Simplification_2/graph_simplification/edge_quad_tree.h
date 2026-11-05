#pragma once

#include <CGAL/Simplification_2/data_structures/quad_tree.h>

#include "utils.h"

namespace CGAL::internal {

	template<class Graph>
	struct EdgeQuadTreeTraits {
		using Element = Graph::Edge_handle;
		using Kernel = Graph::Kernel;

		static Rectangle<Kernel> get_bounding_box(Element elt) {
			Segment<Kernel> seg = elt->curve();
			return utils::boxOf({ seg.start(), seg.end() });
		}

		static bool element_overlaps_rectangle(Element elt, const Rectangle<Kernel>& rect) {
			Segment<Kernel> seg = elt->curve();
			return utils::overlaps(rect, seg);
		}
	};

	template<class Graph>
	using EdgeQuadTree = CGAL::internal::QuadTree<EdgeQuadTreeTraits<Graph>>;
}