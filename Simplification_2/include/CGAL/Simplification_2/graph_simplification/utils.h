#pragma once

#include <CGAL/Simplification_2/core.h>
#include "precision_safe_tests.h"

namespace CGAL::internal::utils {


	template <typename K>
	Rectangle<K> boxOf(const Polygon<K>& p) {
		auto v = p.vertices_begin();
		Number<K> left = v->x();
		Number<K> right = v->x();
		Number<K> bottom = v->y();
		Number<K> top = v->y();

		++v;
		while (v != p.vertices_end()) {
			left = CGAL::min(left, v->x());
			right = CGAL::max(right, v->x());
			bottom = CGAL::min(bottom, v->y());
			top = CGAL::max(top, v->y());
			++v;
		}

		return Rectangle<K>(left, bottom, right, top);
	}

	template <typename K>
	Rectangle<K> boxOf(const std::vector<Point<K>>& p) {
		auto v = p.begin();
		Number<K> left = v->x();
		Number<K> right = v->x();
		Number<K> bottom = v->y();
		Number<K> top = v->y();

		++v;
		while (v != p.end()) {
			left = CGAL::min(left, v->x());
			right = CGAL::max(right, v->x());
			bottom = CGAL::min(bottom, v->y());
			top = CGAL::max(top, v->y());
			++v;
		}

		return Rectangle<K>(left, bottom, right, top);
	}

	template <typename K>
	Rectangle<K> boxOf(const Rectangle<K>& a, const Rectangle<K>& b) {
		Number<K> left = CGAL::min(a.xmin(), b.xmin());
		Number<K> right = CGAL::max(a.xmax(), b.xmax());
		Number<K> bottom = CGAL::min(a.ymin(), b.ymin());
		Number<K> top = CGAL::max(a.ymax(), b.ymax());

		return Rectangle<K>(left, bottom, right, top);
	}

	template <typename K>
	Rectangle<K> boxOf(const Point<K>& a, const Point<K>& b, const Point<K>& c) {
		Number<K> left = CGAL::min(a.x(), CGAL::min(b.x(), c.x()));
		Number<K> right = CGAL::max(a.x(), CGAL::max(b.x(), c.x()));
		Number<K> bottom = CGAL::min(a.y(), CGAL::min(b.y(), c.y()));
		Number<K> top = CGAL::max(a.y(), CGAL::max(b.y(), c.y()));

		return Rectangle<K>(left, bottom, right, top);
	}

	template<typename K>
	Rectangle<K> boxOf(const Triangle<K>& T1, const Triangle<K>& T2) {

		Number<K> left = CGAL::min(CGAL::min(T1[0].x(), CGAL::min(T1[1].x(), T1[2].x())),
			CGAL::min(T2[0].x(), CGAL::min(T2[1].x(), T2[2].x())));
		Number<K> right = CGAL::max(CGAL::max(T1[0].x(), CGAL::max(T1[1].x(), T1[2].x())),
			CGAL::max(T2[0].x(), CGAL::max(T2[1].x(), T2[2].x())));

		Number<K> bottom = CGAL::min(CGAL::min(T1[0].y(), CGAL::min(T1[1].y(), T1[2].y())),
			CGAL::min(T2[0].y(), CGAL::min(T2[1].y(), T2[2].y())));
		Number<K> top = CGAL::max(CGAL::max(T1[0].y(), CGAL::max(T1[1].y(), T1[2].y())),
			CGAL::max(T2[0].y(), CGAL::max(T2[1].y(), T2[2].y())));

		return Rectangle<K>(left, bottom, right, top);
	}

	template<typename K>
	Rectangle<K> boxOf(const std::initializer_list<Point<K>> pts) {

		Number<K> left = 0, right = 0, bottom = 0, top = 0;

		bool first = true;
		for (Point<K> pt : pts) {
			if (first) {
				left = right = pt.x();
				top = bottom = pt.y();
				first = false;
			}
			else {
				if (pt.x() < left) {
					left = pt.x();
				}
				else if (pt.x() > right) {
					right = pt.x();
				}

				if (pt.y() < bottom) {
					bottom = pt.y();
				}
				else if (pt.y() > top) {
					top = pt.y();
				}
			}
		}

		Rectangle<K> box(left, bottom, right, top);
		return box;
	}

	template<class P, typename K>
	Rectangle<K> boxOfWithGetPoint(const std::vector<P>& elements) {

		Number<K> left = 0, right = 0, bottom = 0, top = 0;

		bool first = true;
		for (P elt : elements) {
			Point<K>& pt = elt->getPoint();
			if (first) {
				left = right = pt.x();
				top = bottom = pt.y();
				first = false;
			}
			else {
				if (pt.x() < left) {
					left = pt.x();
				}
				else if (pt.x() > right) {
					right = pt.x();
				}

				if (pt.y() < bottom) {
					bottom = pt.y();
				}
				else if (pt.y() > top) {
					top = pt.y();
				}
			}
		}

		Rectangle<K> box(left, bottom, right, top);
		return box;
	}

	template<typename K>
	bool encloses(const Rectangle<K>& larger, const Rectangle<K>& smaller) {
		return larger.xmin() <= smaller.xmin() && larger.ymin() <= smaller.ymin() &&
			larger.xmax() >= smaller.xmax() && larger.ymax() >= smaller.ymax();
	}

	template<typename K>
	bool disjoint(const Rectangle<K>& a, const Rectangle<K>& b) {
		return a.xmax() < b.xmin() || a.xmin() > b.xmax() || a.ymax() < b.ymin() ||
			a.ymin() > b.ymax();
	}

	template<typename K>
	bool contains(const Rectangle<K>& a, const Point<K>& pt, const Number<K> prec = 0) {
		return a.xmin() - prec <= pt.x() && pt.x() <= a.xmax() + prec && a.ymin() - prec <= pt.y() && pt.y() <= a.ymax() + prec;
	}

	template<typename K>
	bool samePoint(const Point<K>& a, const Point<K>& b, const Number<K> prec = 0) {
		return a.x() - prec <= b.x() && b.x() <= a.x() + prec
			&& a.y() - prec <= b.y() && b.y() <= a.y() + prec;
	}

	template<typename K>
	bool overlaps(const Rectangle<K>& a, const Segment<K>& seg) {
		if (CGAL::intersection(a, seg)) {
			return true;
		}
		else {
			return false;
		}
	}

	template<typename T>
	bool listRemove(T* elt, std::vector<T*>& vec) {

		auto pos = std::find(vec.begin(), vec.end(), elt);
		if (pos != vec.end()) {
			vec.erase(pos);
			return true;
		}
		else {
			return false;
		}
	}

	template<typename T>
	T* swapRemove(int index, std::vector<T*>& vec) {
		if (index == vec.size() - 1) {
			vec.pop_back();
			return nullptr;
		}
		else {
			vec[index] = vec[vec.size() - 1];
			vec.pop_back();
			return vec[index];
		}
	}

	template<typename T>
	void listReplace(T* oldelt, T* newelt, std::vector<T*>& vec) {
		for (int i = 0; i < vec.size(); i++) {
			if (vec[i] == oldelt) {
				vec[i] = newelt;
				return;
			}

		}
		assert(false);
	}

	inline Number<Inexact> solveQuadraticEquationForSmallestPositive(Number<Inexact> a, Number<Inexact> b, Number<Inexact> c) {
		if (safe_test::close(a, 0)) {
			// b x + c = 0 --> x = -c/b
			return -c / b;
		}
		else {
			Number<Inexact> d = b * b - 4 * a * c;
			if (d < -M_EPSILON) {
				return -1;
			}
			else if (d < M_EPSILON) {
				return -b / (2 * a);
			}
			else {
				Number<Inexact> s1 = (-b + std::sqrt(d)) / (2 * a);
				Number<Inexact> s2 = (-b - std::sqrt(d)) / (2 * a);
				if (s1 <= 0) {
					return s2;
				}
				else if (s2 <= 0) {
					return s1;
				}
				else {
					return CGAL::min(s1, s2);
				}
			}
		}
	}
}