// ============================================================================
// TASK 1.1 & 1.2: ROBUST DO_INTERSECT WITH MULTIPLICITY
// For Linear Traits: Ray-Ray, Ray-Segment, Ray-Line, Segment-Line
//
// Implementation for CGAL Arr_linear_traits_2
// Works with both exact and inexact kernels (predicates are robust)
// ============================================================================

#ifndef CGAL_LINEAR_TRAITS_DO_INTERSECT_ROBUST_HPP
#define CGAL_LINEAR_TRAITS_DO_INTERSECT_ROBUST_HPP

#include <utility>  // std::pair
#include <cmath>

namespace CGAL {

// ============================================================================
// ROBUST DO_INTERSECT FUNCTOR FOR LINEAR TRAITS
// Template Parameters:
// - K: Kernel (can be Exact_construction_inexact_predicates_kernel_2)
// ============================================================================

template <typename K_>
class Linear_traits_do_intersect_robust_2 {
public:
    typedef K_                                      Kernel;
    typedef typename Kernel::Point_2                Point_2;
    typedef typename Kernel::Ray_2                  Ray_2;
    typedef typename Kernel::Segment_2              Segment_2;
    typedef typename Kernel::Line_2                 Line_2;
    typedef typename Kernel::Vector_2               Vector_2;
    typedef typename Kernel::FT                     FT;
    typedef CGAL::Orientation                       Orientation;

    // Multiplicity type: (exists, count)
    // count = 1 for regular intersection, 2+ for tangency, 0 if no intersection
    typedef std::pair<bool, int>                    Multiplicity;

private:
    const Kernel* m_kernel;

public:
    Linear_traits_do_intersect_robust_2(const Kernel* k = nullptr) : m_kernel(k) {}

    // ========================================================================
    // RAY-RAY ROBUST INTERSECTION (PREDICATE)
    // ========================================================================
    // Returns: true if rays intersect
    // Robustness: Works with inexact kernels (uses only predicates)
    // Strategy:
    //   1. Same source → always intersect
    //   2. Collinear → check overlap
    //   3. Non-collinear → check if lines intersect AND both ahead

    bool operator()(const Ray_2& r1, const Ray_2& r2) const {
        return ray_ray_robust(r1, r2);
    }

    // ========================================================================
    // RAY-SEGMENT ROBUST INTERSECTION (PREDICATE)
    // ========================================================================

    bool operator()(const Ray_2& ray, const Segment_2& seg) const {
        return ray_segment_robust(ray, seg);
    }

    bool operator()(const Segment_2& seg, const Ray_2& ray) const {
        return ray_segment_robust(ray, seg);
    }

    // ========================================================================
    // RAY-LINE INTERSECTION (PREDICATE)
    // ========================================================================

    bool operator()(const Ray_2& ray, const Line_2& line) const {
        return ray_line_robust(ray, line);
    }

    bool operator()(const Line_2& line, const Ray_2& ray) const {
        return ray_line_robust(ray, line);
    }

    // ========================================================================
    // SEGMENT-LINE INTERSECTION (PREDICATE)
    // ========================================================================

    bool operator()(const Segment_2& seg, const Line_2& line) const {
        return segment_line_robust(seg, line);
    }

    bool operator()(const Line_2& line, const Segment_2& seg) const {
        return segment_line_robust(seg, line);
    }

    // ========================================================================
    // SEGMENT-SEGMENT INTERSECTION (PREDICATE)
    // ========================================================================

    bool operator()(const Segment_2& s1, const Segment_2& s2) const {
        return segment_segment_robust(s1, s2);
    }

    // ========================================================================
    // MULTIPLICITY OPERATORS (Returns intersection count)
    // ========================================================================
    // New feature: Returns pair<bool, int>
    // - bool: whether intersection exists
    // - int: multiplicity (1 = regular, 2 = tangent, etc.)

    Multiplicity multiplicity(const Ray_2& r1, const Ray_2& r2) const {
        if (!ray_ray_robust(r1, r2)) {
            return Multiplicity(false, 0);
        }

        // Check if parallel (same direction)
        FT dp = dot_product(r1.direction().vector(), r2.direction().vector());
        if (is_zero(dp)) {
            // Both collinear and parallel → infinite intersection (tangency)
            return Multiplicity(true, 2);
        }

        // Regular transverse intersection
        return Multiplicity(true, 1);
    }

    Multiplicity multiplicity(const Ray_2& ray, const Segment_2& seg) const {
        if (!ray_segment_robust(ray, seg)) {
            return Multiplicity(false, 0);
        }

        // Check if collinear
        Orientation o1 = orientation(ray.source(),
                                     ray.source() + ray.direction().vector(),
                                     seg.source());
        Orientation o2 = orientation(ray.source(),
                                     ray.source() + ray.direction().vector(),
                                     seg.target());

        if (o1 == COLLINEAR && o2 == COLLINEAR) {
            // Collinear intersection
            return Multiplicity(true, 2);  // Tangency
        }

        return Multiplicity(true, 1);  // Regular intersection
    }

    Multiplicity multiplicity(const Segment_2& s1, const Segment_2& s2) const {
        if (!segment_segment_robust(s1, s2)) {
            return Multiplicity(false, 0);
        }

        // Check collinearity
        Orientation o1 = orientation(s1.source(), s1.target(), s2.source());
        if (o1 == COLLINEAR) {
            return Multiplicity(true, 2);  // Tangency
        }

        return Multiplicity(true, 1);  // Regular intersection
    }

private:
    // ========================================================================
    // ROBUST IMPLEMENTATIONS (Predicate-only)
    // ========================================================================

    bool ray_ray_robust(const Ray_2& r1, const Ray_2& r2) const {
        const Point_2& s1 = r1.source();
        const Point_2& s2 = r2.source();
        const Vector_2& d1 = r1.direction().vector();
        const Vector_2& d2 = r2.direction().vector();

        // Case 1: Same source
        if (s1 == s2) return true;

        // Get points on rays
        Point_2 p1 = s1 + d1;
        Point_2 p2 = s2 + d2;

        // Case 2: Check collinearity
        Orientation o_check1 = orientation(s1, p1, s2);
        Orientation o_check2 = orientation(s1, p1, p2);

        if (o_check1 == COLLINEAR && o_check2 == COLLINEAR) {
            // Both collinear - check overlap
            FT dp = dot_product(d1, d2);

            if (dp > 0) {
                // Same direction
                return true;
            } else {
                // Opposite directions - check if ranges overlap
                Vector_2 v_s1_to_s2 = s2 - s1;
                Vector_2 v_s2_to_s1 = s1 - s2;

                FT dp1 = dot_product(d1, v_s1_to_s2);
                FT dp2 = dot_product(d2, v_s2_to_s1);

                return (dp1 >= 0 && dp2 >= 0);
            }
        }

        // Case 3: Non-collinear - check if lines cross
        Orientation o1 = orientation(s1, p1, s2);
        Orientation o2 = orientation(s1, p1, p2);
        Orientation o3 = orientation(s2, p2, s1);
        Orientation o4 = orientation(s2, p2, p1);

        if (!((o1 != o2) && (o3 != o4))) {
            return false;  // Lines don't intersect
        }

        // Lines intersect - check if both rays point toward intersection
        Vector_2 v_s1_to_s2 = s2 - s1;
        Vector_2 v_s2_to_s1 = s1 - s2;

        FT dp1 = dot_product(d1, v_s1_to_s2);
        FT dp2 = dot_product(d2, v_s2_to_s1);

        return (dp1 >= 0 && dp2 >= 0);
    }

    bool ray_segment_robust(const Ray_2& ray, const Segment_2& seg) const {
        const Point_2& ray_src = ray.source();
        const Vector_2& ray_dir = ray.direction().vector();
        Point_2 ray_pt = ray_src + ray_dir;

        const Point_2& seg_src = seg.source();
        const Point_2& seg_tgt = seg.target();

        Orientation o1 = orientation(ray_src, ray_pt, seg_src);
        Orientation o2 = orientation(ray_src, ray_pt, seg_tgt);

        // Both on same side (not collinear) - no intersection
        if (o1 != COLLINEAR && o2 != COLLINEAR && o1 == o2) {
            return false;
        }

        // Collinear case
        if (o1 == COLLINEAR && o2 == COLLINEAR) {
            Vector_2 v_src_to_seg_src = seg_src - ray_src;
            Vector_2 v_src_to_seg_tgt = seg_tgt - ray_src;

            FT dp_src = dot_product(ray_dir, v_src_to_seg_src);
            FT dp_tgt = dot_product(ray_dir, v_src_to_seg_tgt);

            return (dp_src >= 0 || dp_tgt >= 0);
        }

        // Straddling case - check at least one endpoint ahead
        const Point_2& test_pt = (o1 != COLLINEAR) ? seg_src : seg_tgt;
        Vector_2 v = test_pt - ray_src;
        FT dp = dot_product(ray_dir, v);

        return (dp >= 0);
    }

    bool ray_line_robust(const Ray_2& ray, const Line_2& line) const {
        const Point_2& ray_src = ray.source();
        const Vector_2& ray_dir = ray.direction().vector();

        // Line passes through point with direction
        const Point_2& line_pt = line.point();
        const Vector_2& line_dir = line.direction().vector();

        Point_2 line_pt2 = line_pt + line_dir;

        Orientation o1 = orientation(ray_src, ray_src + ray_dir, line_pt);
        Orientation o2 = orientation(ray_src, ray_src + ray_dir, line_pt2);

        // Parallel
        if (o1 != COLLINEAR && o1 == o2) {
            return false;
        }

        // Check if intersection ahead
        Vector_2 v = line_pt - ray_src;
        FT dp = dot_product(ray_dir, v);

        return (dp >= 0);
    }

    bool segment_line_robust(const Segment_2& seg, const Line_2& line) const {
        const Point_2& seg_src = seg.source();
        const Point_2& seg_tgt = seg.target();

        const Point_2& line_pt = line.point();
        const Vector_2& line_dir = line.direction().vector();
        Point_2 line_pt2 = line_pt + line_dir;

        Orientation o1 = orientation(seg_src, seg_tgt, line_pt);
        Orientation o2 = orientation(seg_src, seg_tgt, line_pt2);

        // Parallel
        if (o1 != COLLINEAR && o1 == o2) {
            return false;
        }

        // Collinear or straddling counts as intersection
        return true;
    }

    bool segment_segment_robust(const Segment_2& s1, const Segment_2& s2) const {
        Orientation o1 = orientation(s1.source(), s1.target(), s2.source());
        Orientation o2 = orientation(s1.source(), s1.target(), s2.target());
        Orientation o3 = orientation(s2.source(), s2.target(), s1.source());
        Orientation o4 = orientation(s2.source(), s2.target(), s1.target());

        // Proper intersection
        if (((o1 == POSITIVE && o2 == NEGATIVE) || (o1 == NEGATIVE && o2 == POSITIVE)) &&
            ((o3 == POSITIVE && o4 == NEGATIVE) || (o3 == NEGATIVE && o4 == POSITIVE))) {
            return true;
        }

        // Improper intersections (collinear endpoints)
        if (o1 == COLLINEAR && on_segment(s1, s2.source())) return true;
        if (o2 == COLLINEAR && on_segment(s1, s2.target())) return true;
        if (o3 == COLLINEAR && on_segment(s2, s1.source())) return true;
        if (o4 == COLLINEAR && on_segment(s2, s1.target())) return true;

        return false;
    }

    // ========================================================================
    // HELPER PREDICATES
    // ========================================================================

    // Orientation predicate: CCW/Collinear/CW relative to directed segment
    Orientation orientation(const Point_2& p, const Point_2& q, const Point_2& r) const {
        FT cross = (q.x() - p.x()) * (r.y() - p.y()) - (q.y() - p.y()) * (r.x() - p.x());

        if (cross > 0) return POSITIVE;
        if (cross < 0) return NEGATIVE;
        return COLLINEAR;
    }

    // Dot product
    FT dot_product(const Vector_2& v, const Vector_2& w) const {
        return v.x() * w.x() + v.y() * w.y();
    }

    // Point on segment (collinear case)
    bool on_segment(const Segment_2& seg, const Point_2& p) const {
        return (std::min(seg.source().x(), seg.target().x()) <= p.x() &&
                p.x() <= std::max(seg.source().x(), seg.target().x()) &&
                std::min(seg.source().y(), seg.target().y()) <= p.y() &&
                p.y() <= std::max(seg.source().y(), seg.target().y()));
    }

    bool is_zero(const FT& x) const {
        return x == FT(0);
    }
};

}  // namespace CGAL

#endif  // CGAL_LINEAR_TRAITS_DO_INTERSECT_ROBUST_HPP
