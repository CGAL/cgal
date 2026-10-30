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

template <PointQuadTreeTraits PQT> struct PQTNode {

	using Node = PQTNode<PQT>;
	using Rect = Rectangle<typename PQT::Kernel>;
	using Element = PQT::Element;

	Node* parent;
	Node* lt;
	Node* lb;
	Node* rt;
	Node* rb;
	const Rect rect;
	std::vector<Element>* elts;

	PQTNode(Node* p, Rect r, const bool is_leaf) : parent(p), rect(r) {
		lt = nullptr;
		lb = nullptr;
		rt = nullptr;
		rb = nullptr;
		elts = is_leaf ? new std::vector<Element>() : nullptr;
	}

	~PQTNode() {
		delete lt;
		delete lb;
		delete rt;
		delete rb;
		delete elts;
	}
};
} // namespace detail

template <PointQuadTreeTraits PQT>
void PointQuadTree<PQT>::findContainedRecursive(Node* n, const Rectangle<Kernel>& query,
                                                ElementCallback act) {

	if (n == nullptr) {
		return;
	} else if (disjoint(n->rect, query)) {
		return;
	}

	if (n->elts == nullptr) {
		// internal node
		findContainedRecursive(n->lb, query, act);
		findContainedRecursive(n->lt, query, act);
		findContainedRecursive(n->rb, query, act);
		findContainedRecursive(n->rt, query, act);
	} else if (encloses(query, n->rect)) {
		for (Element elt : *(n->elts)) {
			act(elt);
		}
	} else {
		for (Element elt : *(n->elts)) {
			if (contains(query, PQT::get_point(elt), 0)) {
				act(elt);
			}
		}
	}
}
template <PointQuadTreeTraits PQT> Rectangle<typename PQT::Kernel> PointQuadTree<PQT>::root_box() {
	return root->rect;
}

template <PointQuadTreeTraits PQT>
template <bool extend>
detail::PQTNode<PQT>* PointQuadTree<PQT>::find(Element elt) {
	int d = 0;
	Node* n = root;
	const Point<Kernel>& pt = PQT::get_point(elt);

	while (n->elts == nullptr) {
		d++;

		Number<Kernel> xmid = (n->rect.xmin() + n->rect.xmax()) / 2;
		Number<Kernel> ymid = (n->rect.ymin() + n->rect.ymax()) / 2;

		if (pt.x() <= xmid) {
			// left
			if (pt.y() <= ymid) {
				// bottom
				if (n->lb == nullptr) {
					if constexpr (extend) {
						n->lb =
						    new Node(n, Rectangle<Kernel>(n->rect.xmin(), n->rect.ymin(), xmid, ymid),
						             d >= maxdepth);
					} else {
						return nullptr;
					}
				}
				n = n->lb;
			} else {
				// top
				if (n->lt == nullptr) {
					if constexpr (extend) {
						n->lt =
						    new Node(n, Rectangle<Kernel>(n->rect.xmin(), ymid, xmid, n->rect.ymax()),
						             d >= maxdepth);
					} else {
						return nullptr;
					}
				}
				n = n->lt;
			}
		} else {
			// right
			if (pt.y() <= ymid) {
				// bottom
				if (n->rb == nullptr) {
					if constexpr (extend) {
						n->rb =
						    new Node(n, Rectangle<Kernel>(xmid, n->rect.ymin(), n->rect.xmax(), ymid),
						             d >= maxdepth);
					} else {
						return nullptr;
					}
				}
				n = n->rb;
			} else {
				// top
				if (n->rt == nullptr) {
					if constexpr (extend) {
						n->rt =
						    new Node(n, Rectangle<Kernel>(xmid, ymid, n->rect.xmax(), n->rect.ymax()),
						             d >= maxdepth);
					} else {
						return nullptr;
					}
				}
				n = n->rt;
			}
		}
	}

	return n;
}

template <PointQuadTreeTraits PQT>
PointQuadTree<PQT>::PointQuadTree(const Rectangle<Kernel>& box, int depth) {
	root = new Node(nullptr, box, depth == 0);
	maxdepth = depth;
}

template <PointQuadTreeTraits PQT> PointQuadTree<PQT>::~PointQuadTree() {
	delete root;
}

template <PointQuadTreeTraits PQT> void PointQuadTree<PQT>::clear() {
	// just delete all nodes except the root
	delete root->lb;
	delete root->rb;
	delete root->lt;
	delete root->rt;
	root->lb = nullptr;
	root->rb = nullptr;
	root->lt = nullptr;
	root->rt = nullptr;
}

template <PointQuadTreeTraits PQT> void PointQuadTree<PQT>::insert(Element elt) {
	Node* n = find<true>(elt);
	n->elts->push_back(elt);
}

template <PointQuadTreeTraits PQT> bool PointQuadTree<PQT>::remove(Element elt) {
	Node* n = find<false>(elt);

	if (n == nullptr) {
		// it's supposed leaf doesn't exist
		return false;
	}

	auto position = std::find(n->elts->begin(), n->elts->end(), elt);
	if (position == n->elts->end()) {
		// the node doesn't actually contain the element
		return false;
	}

	// it's contained: erase and clean up the tree if possible
	n->elts->erase(position);

	// clean up the tree if possible
	if (n->elts->empty()) {
		do {
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
		} while (n != root && n->lb == nullptr && n->lt == nullptr && n->rb == nullptr &&
		         n->rt == nullptr);
	}

	return true;
}

template <PointQuadTreeTraits PQT>
void PointQuadTree<PQT>::findContained(const Rectangle<Kernel>& query, ElementCallback act) {
	findContainedRecursive(root, query, act);
}

template <PointQuadTreeTraits PQT>
PQT::Element PointQuadTree<PQT>::findElementRecursive(Node* n, const Point<Kernel>& query,
                                                      const Number<Kernel> prec) {
	if (n == nullptr || !contains(n->rect, query, prec)) {
		// outside of the node's rectangle
		return nullptr;
	} else if (n->elts == nullptr) {
		// internal node
		Element elt = findElementRecursive(n->lb, query, prec);
		if (elt != nullptr) {
			return elt;
		}
		elt = findElementRecursive(n->lt, query, prec);
		if (elt != nullptr) {
			return elt;
		}
		elt = findElementRecursive(n->rb, query, prec);
		if (elt != nullptr) {
			return elt;
		}
		elt = findElementRecursive(n->rt, query, prec);
		return elt;
	} else {
		// leaf node
		for (Element elt : *(n->elts)) {
			if (same_point(query, PQT::get_point(elt), prec)) {
				return elt;
			}
		}
		return nullptr;
	}
}

template <PointQuadTreeTraits PQT>
PQT::Element PointQuadTree<PQT>::findElement(const Point<Kernel>& query, const Number<Kernel> prec) {
	return findElementRecursive(root, query, prec);
}

template <PointQuadTreeTraits PQT>
bool PointQuadTree<PQT>::encloses(const Rectangle<Kernel>& larger, const Rectangle<Kernel>& smaller) {
	return larger.xmin() <= smaller.xmin() && larger.ymin() <= smaller.ymin() &&
	       larger.xmax() >= smaller.xmax() && larger.ymax() >= smaller.ymax();
}

template <PointQuadTreeTraits PQT>
bool PointQuadTree<PQT>::disjoint(const Rectangle<Kernel>& a, const Rectangle<Kernel>& b) {
	return a.xmax() < b.xmin() || a.xmin() > b.xmax() || a.ymax() < b.ymin() || a.ymin() > b.ymax();
}

template <PointQuadTreeTraits PQT>
bool PointQuadTree<PQT>::contains(const Rectangle<Kernel>& rect, const Point<Kernel>& point,
                                  const Number<Kernel> prec) {
	return rect.xmin() - prec <= point.x() && point.x() <= rect.xmax() + prec &&
	       rect.ymin() - prec <= point.y() && point.y() <= rect.ymax() + prec;
}

template <PointQuadTreeTraits PQT>
bool PointQuadTree<PQT>::same_point(const Point<Kernel>& a, const Point<Kernel>& b,
                                    const Number<Kernel> prec) {
	return a.x() - prec <= b.x() && b.x() <= a.x() + prec && a.y() - prec <= b.y() &&
	       b.y() <= a.y() + prec;
}

} // namespace CGAL::internal