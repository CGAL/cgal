#pragma once

#include <CGAL/Simplification_2/core.h>

namespace CGAL::internal::safe_test {

	inline bool close(const Number<Exact> a, const Number<Exact> b) {
		return a == b;
	}
	inline bool close(const Number<Inexact> a, const Number<Inexact> b) {
		return std::abs(a - b) < M_EPSILON;
	}

	inline bool leq(const Number<Exact> a, const Number<Exact> b) {
		return a <= b;
	}
	inline bool leq(const Number<Inexact> a, const Number<Inexact> b) {
		return a <= b + M_EPSILON;
	}

	inline bool same_point(const Point<Exact> a, const Point<Exact> b) {
		return a == b;
	}
	inline bool same_point(const Point<Inexact> a, const Point<Inexact> b) {
		return close(a.x(), b.x()) && close(a.y(), b.y());
	}

	inline bool aligned_vectors(const Vector<Exact> a, const Vector<Exact> b) {
		return a.direction() == b.direction();
	}
	inline bool aligned_vectors(const Vector<Inexact> a, const Vector<Inexact> b) {
		Vector<Inexact> a_norm = a / std::sqrt(a.squared_length());
		Vector<Inexact> b_norm = b / std::sqrt(b.squared_length());
		return close(CGAL::determinant(a_norm, b_norm), 0) && CGAL::scalar_product(a, b) > 0;
	}

	inline bool aligned_or_opposite_vectors(const Vector<Exact> a, const Vector<Exact> b) {
		return a.direction() == b.direction() || a.direction() == -b.direction();
	}
	inline bool aligned_or_opposite_vectors(const Vector<Inexact> a, const Vector<Inexact> b) {
		Vector<Inexact> a_norm = a / std::sqrt(a.squared_length());
		Vector<Inexact> b_norm = b / std::sqrt(b.squared_length());
		return close(CGAL::determinant(a_norm, b_norm), 0);
	}

	inline bool collinear(const Point<Exact> a, const Point<Exact> b, const Point<Exact> c) {
		return CGAL::collinear(a, b, c);
	}
	inline bool collinear(const Point<Inexact> a, const Point<Inexact> b, const Point<Inexact> c) {
		return aligned_or_opposite_vectors(b - a, c - a);
	}

	inline bool point_on_line(const Point<Exact> p, const Line<Exact> l) {
		return l.has_on_boundary(p);
	}
	inline bool point_on_line(const Point<Inexact> p, const Line<Inexact> l) {
		return CGAL::squared_distance(l, p) < M_EPSILON;
	}

	inline bool convex_contains(const Point<Exact> a, const Polygon<Exact> p) {
		return !p.has_on_unbounded_side(a);
	}
	inline bool convex_contains(const Point<Inexact> a, const Polygon<Inexact> p) {
		const size_t n = p.size();

		Vector<Inexact> dirPrev = p[n - 1] - a;
		dirPrev = dirPrev / std::sqrt(dirPrev.squared_length());
		Vector<Inexact> dirCurr = p[0] - a;
		dirCurr = dirCurr / std::sqrt(dirCurr.squared_length());

		Number<Inexact> cp = CGAL::determinant(dirPrev, dirCurr);
		auto sig = CGAL::sign(cp);

		if (close(cp, 0)) {
			// either vertex is on line segment (inside) or in parallel with it (outside)
			// inside == dotproduct <= 0
			return leq(CGAL::scalar_product(dirPrev, dirCurr), 0);
		}
		else {
			for (size_t i = 1; i < n; i++) {
				dirPrev = dirCurr;
				dirCurr = p[i] - a;
				dirCurr = dirCurr / std::sqrt(dirCurr.squared_length());

				cp = CGAL::determinant(dirPrev, dirCurr);
				auto s = CGAL::sign(cp);
				if (close(cp, 0)) {
					// either vertex is on line segment (inside) or in parallel with it (outside)
					// inside == dotproduct <= 0
					return leq(CGAL::scalar_product(dirPrev, dirCurr), 0);
				}
				else if (sig == s) {
					// continue, same sign
				}
				else {
					// opposite signs
					return false;
				}
			}

			return true;
		}
	}

}

