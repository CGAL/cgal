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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF TEMPLATE FUNCTIONS
// Do not include this file, but the .h file instead
// -----------------------------------------------------------------------------

namespace CGAL::internal {

namespace detail {
template <QuadTreeTraits QT> struct QTNode {

	using Node = QTNode<QT>;
	using Rect = Rectangle<typename QT::Kernel>;
	using Element = QT::Element;

	Node* parent;
	Node* lt = nullptr;
	Node* lb = nullptr;
	Node* rt = nullptr;
	Node* rb = nullptr;
	Rect rect;
	bool inf_left, inf_right, inf_bottom, inf_top;
	std::vector<Element> elts;

	QTNode(Node* parent, Rect rect, bool inf_left, bool inf_right, bool inf_bottom, bool inf_top)
	    : parent(parent), rect(rect), inf_left(inf_left), inf_right(inf_right),
	      inf_bottom(inf_bottom), inf_top(inf_top) {}

	~QTNode() {
		delete lt;
		delete lb;
		delete rt;
		delete rb;
	}
};
} // namespace detail

template <QuadTreeTraits QT> Rectangle<typename QT::Kernel> QuadTree<QT>::root_box() {
	return root->rect;
}

template <QuadTreeTraits QT>
bool QuadTree<QT>::encloses(const Rectangle<Kernel>& rect, Node* node) {
	// node extends further to the left
	if (node->inf_left || rect.xmin() > node->rect.xmin()) {
		return false;
	}

	// node extends further to the right
	if (node->inf_right || rect.xmax() < node->rect.xmax()) {
		return false;
	}

	// node extends further to the bottom
	if (node->inf_bottom || rect.ymin() > node->rect.ymin()) {
		return false;
	}

	// node extends further to the top
	if (node->inf_top || rect.ymax() < node->rect.ymax()) {
		return false;
	}

	return true;
}

template <QuadTreeTraits QT>
bool QuadTree<QT>::disjoint(Node* node, const Rectangle<Kernel>& rect) {
	// rect left of node
	if (!node->inf_left && rect.xmax() < node->rect.xmin()) {
		return true;
	}

	// rect right of node
	if (!node->inf_right && rect.xmin() > node->rect.xmax()) {
		return true;
	}

	// rect below node
	if (!node->inf_bottom && rect.ymax() < node->rect.ymin()) {
		return true;
	}

	// rect above node
	if (!node->inf_top && rect.ymin() > node->rect.ymax()) {
		return true;
	}

	return false;
}

template <QuadTreeTraits QT>
void QuadTree<QT>::findOverlappedRecursive(Node* n, const Rectangle<Kernel>& query,
                                           ElementCallback act) {

	if (n == nullptr) {
		return;
	} else if (disjoint(n, query)) {
		return;
	} else if (encloses(query, n)) {
		for (Element elt : n->elts) {
			act(elt);
		}
	} else {
		for (Element elt : n->elts) {
			if (QT::element_overlaps_rectangle(elt, query)) {
				act(elt);
			}
		}
	}

	findOverlappedRecursive(n->lb, query, act);
	findOverlappedRecursive(n->lt, query, act);
	findOverlappedRecursive(n->rb, query, act);
	findOverlappedRecursive(n->rt, query, act);
}

template <QuadTreeTraits QT>
template <bool extend>
detail::QTNode<QT>* QuadTree<QT>::find(Element elt) {
	int d = 0;
	Node* n = root;

	Rectangle<Kernel> bb = QT::get_bounding_box(elt);

	// NB: extending, n will never become null
	while (d < maxdepth) {
		d++;

		Number<Kernel> xmid = (n->rect.xmin() + n->rect.xmax()) / 2;
		Number<Kernel> ymid = (n->rect.ymin() + n->rect.ymax()) / 2;
		Number<Kernel> fx = (n->rect.xmax() - n->rect.xmin()) / 2.0 * fuzziness;
		Number<Kernel> fy = (n->rect.ymax() - n->rect.ymin()) / 2.0 * fuzziness;

		if (bb.xmax() <= xmid + fx) {
			// left
			if (bb.ymax() <= ymid + fy) {
				// bottom
				if (n->lb == nullptr) {
					if constexpr (extend) {
						n->lb = new Node(
						    n,
						    Rectangle<Kernel>(n->rect.xmin(), n->rect.ymin(), xmid + fx, ymid + fy),
						    n->inf_left, false, n->inf_bottom, false);
					} else {
						return nullptr;
					}
				}
				n = n->lb;
			} else if (bb.ymin() >= ymid - fy) {
				// top
				if (n->lt == nullptr) {
					if constexpr (extend) {
						n->lt = new Node(
						    n,
						    Rectangle<Kernel>(n->rect.xmin(), ymid - fy, xmid + fx, n->rect.ymax()),
						    n->inf_left, false, false, n->inf_top);
					} else {
						return nullptr;
					}
				}
				n = n->lt;
			} else {
				// here
				return n;
			}
		} else if (bb.xmin() >= xmid - fx) {
			// right
			if (bb.ymax() <= ymid + fy) {
				// bottom
				if (n->rb == nullptr) {
					if constexpr (extend) {
						n->rb = new Node(
						    n,
						    Rectangle<Kernel>(xmid - fx, n->rect.ymin(), n->rect.xmax(), ymid + fy),
						    false, n->inf_right, n->inf_bottom, false);
					} else {
						return nullptr;
					}
				}
				n = n->rb;
			} else if (bb.ymin() >= ymid - fy) {
				// top
				if (n->rt == nullptr) {
					if constexpr (extend) {
						n->rt = new Node(
						    n,
						    Rectangle<Kernel>(xmid - fx, ymid - fy, n->rect.xmax(), n->rect.ymax()),
						    false, n->inf_right, false, n->inf_top);
					} else {
						return nullptr;
					}
				}
				n = n->rt;
			} else {
				// here
				return n;
			}
		} else {
			// here
			return n;
		}
	}

	return n;
}

template <QuadTreeTraits QT>
QuadTree<QT>::QuadTree(const Rectangle<Kernel>& box, int depth, Number<Kernel> fuzz) {
	root = new Node(nullptr, box, true, true, true, true);
	maxdepth = depth;
	fuzziness = fuzz;
}

template <QuadTreeTraits QT> QuadTree<QT>::~QuadTree() {
	delete root;
}

template <QuadTreeTraits QT> void QuadTree<QT>::clear() {
	// just delete all nodes except the root
	delete root->lb;
	delete root->rb;
	delete root->lt;
	delete root->rt;
	root->lb = nullptr;
	root->rb = nullptr;
	root->lt = nullptr;
	root->rt = nullptr;

	// and empty the elements at the root
	root->elts.clear();
}

template <QuadTreeTraits QT> void QuadTree<QT>::insert(Element elt) {
	Node* n = find<true>(elt);
	n->elts.push_back(elt);
}

template <QuadTreeTraits QT> bool QuadTree<QT>::remove(Element elt) {
	Node* n = find<false>(elt);

	if (n == nullptr) {
		// it's supposed containing node doesn't exist
		return false;
	}

	auto position = std::find(n->elts.begin(), n->elts.end(), elt);
	if (position == n->elts.end()) {
		// the node doesn't actually contain the element
		return false;
	}

	// it's contained: erase and clean up the tree if possible
	n->elts.erase(position);

	while (n != root && n->elts.empty() && n->lb == nullptr && n->lt == nullptr &&
	       n->rb == nullptr && n->rt == nullptr) {
		Node* p = n->parent;
		if (p->lb == n) {
			p->lb = nullptr;
		} else if (p->lt == n) {
			p->lt = nullptr;
		} else if (p->rb == n) {
			p->rb = nullptr;
		} else if (p->rt == n) {
			p->rt = nullptr;
		}
		delete n;
		n = p;
	}

	return true;
}

template <QuadTreeTraits QT>
void QuadTree<QT>::findOverlapped(const Rectangle<Kernel>& query, ElementCallback act) {
	findOverlappedRecursive(root, query, act);
}

} // namespace CGAL::internal