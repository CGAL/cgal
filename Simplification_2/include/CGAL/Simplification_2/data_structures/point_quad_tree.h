/*
Copyright (C) 2026  TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include <CGAL/Simplification_2/core.h>

namespace CGAL::internal {

template <typename PQT> concept PointQuadTreeTraits = requires(typename PQT::Element elt) {

	typename PQT::Element;
	typename PQT::Kernel;

	{PQT::get_point(elt)}->std::convertible_to<const Point<typename PQT::Kernel>&>;
};

namespace detail {
template <PointQuadTreeTraits PQT> struct PQTNode;
}

template <PointQuadTreeTraits PQT> class PointQuadTree {
  public:
	using Kernel = PQT::Kernel;
	using Element = PQT::Element;
	using ElementCallback = std::function<void(Element)>;

	PointQuadTree(const Rectangle<Kernel>& box, int depth);
	~PointQuadTree();

	void clear();
	void insert(Element elt);
	bool remove(Element elt);

	void findContained(const Rectangle<Kernel>& query, ElementCallback act);
	Element findElement(const Point<Kernel>& query, const Number<Kernel> prec = 0);

    Rectangle<Kernel> root_box();

  private:
	using Node = detail::PQTNode<PQT>;

	Node* root;
	int maxdepth;

	Element findElementRecursive(Node* n, const Point<Kernel>& query, const Number<Kernel> prec);
	void findContainedRecursive(Node* n, const Rectangle<Kernel>& query, ElementCallback act);
	template <bool extend> Node* find(Element elt);

	// Does the larger rectangle enclose the smaller?
	static bool encloses(const Rectangle<Kernel>& larger, const Rectangle<Kernel>& smaller);
	// Are these rectangles disjoint?
	static bool disjoint(const Rectangle<Kernel>& a, const Rectangle<Kernel>& b);
	// Does the rectangle contain the point?
	static bool contains(const Rectangle<Kernel>& rect, const Point<Kernel>& point,
	                     const Number<Kernel> prec);
	// Are these points the same, up to the given precision?
	static bool same_point(const Point<Kernel>& a, const Point<Kernel>& b, const Number<Kernel> prec);
};

} // namespace CGAL::internal
#include "point_quad_tree.hpp"