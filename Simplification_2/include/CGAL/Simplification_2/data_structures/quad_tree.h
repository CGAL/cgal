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

	template <typename QT> concept QuadTreeTraits = requires(typename QT::Element elt, const Rectangle<typename QT::Kernel>& rect) {

		typename QT::Element;
		typename QT::Kernel;

		{
			QT::get_bounding_box(elt)
		} -> std::convertible_to<Rectangle<typename QT::Kernel>>;

		{
			QT::element_overlaps_rectangle(elt,rect)
		} -> std::convertible_to<bool>;
	};

	namespace detail {
		template <QuadTreeTraits QT> 
		class QTNode;
	}

	template <QuadTreeTraits QT> class QuadTree {
	public:
		using Kernel = QT::Kernel;
		using Element = QT::Element;
		using ElementCallback = std::function<void(Element)>;

		QuadTree(const Rectangle<Kernel>& box, int depth, Number<Kernel> fuzz);
		~QuadTree();

		void clear();
		void insert(Element elt);
		bool remove(Element elt);

		void findOverlapped(const Rectangle<Kernel>& query, ElementCallback act);

	    Rectangle<Kernel> root_box();

	private:
		using Node = detail::QTNode<QT>;

		Node* root;
		int maxdepth;
		Number<Kernel> fuzziness;

		// Does the rectangle enclose the (possibly infinite) node?
	    bool encloses(const Rectangle<Kernel>& rect, Node* node);

		// Is the (possibly infinite) node disjoint from the rectangle
	    bool disjoint(Node* node, const Rectangle<Kernel>& rect);

		void findOverlappedRecursive(Node* n, const Rectangle<Kernel>& query, ElementCallback act);

		template <bool extend> Node* find(Element elt);
	};

} // namespace CGAL::internal

#include "quad_tree.hpp"