// ============================================================================
// TASK 2: CIRCLE-SEGMENT TRAITS EXTENSION FOR RAYS AND LINES
// Add support for AosOpenBoundaryTraits_2
//
// This extends the existing circle-segment traits to support:
// - CGAL::Ray_2 in addition to segments
// - CGAL::Line_2 in addition to segments
// ============================================================================

#ifndef CGAL_CIRCLE_SEGMENT_RAY_LINE_EXTENSION_HPP
#define CGAL_CIRCLE_SEGMENT_RAY_LINE_EXTENSION_HPP

namespace CGAL {

// ============================================================================
// OPEN CURVE TYPES (Rays and Lines)
// ============================================================================

template <typename Kernel_>
class Circle_ray_2 {
public:
    typedef Kernel_ Kernel;
    typedef typename Kernel::Circle_2 Circle_2;
    typedef typename Kernel::Ray_2 Ray_2;
    typedef typename Kernel::Point_2 Point_2;

private:
    Circle_2 m_circle;
    Ray_2 m_ray;

public:
    Circle_ray_2() {}
    Circle_ray_2(const Circle_2& c, const Ray_2& r) : m_circle(c), m_ray(r) {}

    const Circle_2& circle() const { return m_circle; }
    const Ray_2& ray() const { return m_ray; }
};

template <typename Kernel_>
class Circle_line_2 {
public:
    typedef Kernel_ Kernel;
    typedef typename Kernel::Circle_2 Circle_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Point_2 Point_2;

private:
    Circle_2 m_circle;
    Line_2 m_line;

public:
    Circle_line_2() {}
    Circle_line_2(const Circle_2& c, const Line_2& l) : m_circle(c), m_line(l) {}

    const Circle_2& circle() const { return m_circle; }
    const Line_2& line() const { return m_line; }
};

// ============================================================================
// DO_INTERSECT FUNCTOR EXTENSION FOR RAYS AND LINES
// ============================================================================

template <typename Kernel_>
class Circle_segment_do_intersect_extended {
public:
    typedef Kernel_ Kernel;
    typedef typename Kernel::Circle_2 Circle_2;
    typedef typename Kernel::Ray_2 Ray_2;
    typedef typename Kernel::Line_2 Line_2;
    typedef typename Kernel::Segment_2 Segment_2;
    typedef typename Kernel::Point_2 Point_2;
    typedef typename Kernel::FT FT;

    // ========================================================================
    // CIRCLE-RAY INTERSECTION
    // ========================================================================
    // Returns: true if circle and ray intersect
    // Strategy: Circle-line intersection + check if on ray

    bool operator()(const Circle_2& circle, const Ray_2& ray) const {
        return circle_ray_intersect(circle, ray);
    }

    bool operator()(const Ray_2& ray, const Circle_2& circle) const {
        return circle_ray_intersect(circle, ray);
    }

    // ========================================================================
    // CIRCLE-LINE INTERSECTION
    // ========================================================================
    // Returns: true if circle and line intersect

    bool operator()(const Circle_2& circle, const Line_2& line) const {
        return circle_line_intersect(circle, line);
    }

    bool operator()(const Line_2& line, const Circle_2& circle) const {
        return circle_line_intersect(circle, line);
    }

    // ========================================================================
    // CIRCLE-SEGMENT (existing, for reference)
    // ========================================================================

    bool operator()(const Circle_2& circle, const Segment_2& seg) const {
        return circle_segment_intersect(circle, seg);
    }

    bool operator()(const Segment_2& seg, const Circle_2& circle) const {
        return circle_segment_intersect(circle, seg);
    }

private:
    // ========================================================================
    // ROBUST CIRCLE-LINE INTERSECTION (core predicate)
    // ========================================================================
    // Uses only distance computation (no root finding)
    // Works with inexact kernels

    bool circle_line_intersect(const Circle_2& circle, const Line_2& line) const {
        // Get line point and direction
        const Point_2& line_pt = line.point();
        const typename Kernel::Direction_2& line_dir = line.direction();

        // Compute distance from circle center to line
        const Point_2& center = circle.center();
        FT radius_sq = circle.squared_radius();

        // Vector from center to line point
        typename Kernel::Vector_2 v = line_pt - center;

        // Perpendicular distance to line: |v - (v·d)d|²
        // where d is unit direction (normalized)
        // For robustness, use (v × d)² / |d|²
        typename Kernel::Vector_2 d = line_dir.vector();

        // Cross product magnitude: |v × d|
        FT cross = v.x() * d.y() - v.y() * d.x();
        FT d_norm_sq = d.x() * d.x() + d.y() * d.y();

        // Distance² from center to line = cross² / d_norm²
        FT dist_sq = (cross * cross) / d_norm_sq;

        // Intersect iff dist² <= radius²
        return (dist_sq <= radius_sq);
    }

    // ========================================================================
    // CIRCLE-RAY INTERSECTION
    // ========================================================================
    // Strategy:
    //   1. Check intersection with infinite line through ray
    //   2. Check if intersection points are ahead of ray source

    bool circle_ray_intersect(const Circle_2& circle, const Ray_2& ray) const {
        // First: do ray and circle's baseline line intersect?
        typename Kernel::Line_2 baseline_line(ray.source(), ray.direction());

        if (!circle_line_intersect(circle, baseline_line)) {
            return false;
        }

        // Second: check if any intersection is ahead of ray source
        const Point_2& center = circle.center();
        const Point_2& ray_src = ray.source();
        const typename Kernel::Vector_2& ray_dir = ray.direction().vector();

        // Vector from ray source to center
        typename Kernel::Vector_2 v = center - ray_src;

        // Projection of (center - source) onto ray direction
        FT projection = v.x() * ray_dir.x() + v.y() * ray_dir.y();

        // If projection >= 0, ray might intersect
        return (projection >= 0);
    }

    // ========================================================================
    // CIRCLE-SEGMENT INTERSECTION (existing)
    // ========================================================================
    // Check if any segment point is within circle radius

    bool circle_segment_intersect(const Circle_2& circle, const Segment_2& seg) const {
        const Point_2& center = circle.center();
        FT radius_sq = circle.squared_radius();

        // Check segment endpoints
        if (distance_squared(center, seg.source()) <= radius_sq) {
            return true;
        }
        if (distance_squared(center, seg.target()) <= radius_sq) {
            return true;
        }

        // Check closest point on segment to circle center
        typename Kernel::Vector_2 seg_vec = seg.target() - seg.source();
        typename Kernel::Vector_2 v = center - seg.source();

        FT seg_len_sq = seg_vec.x() * seg_vec.x() + seg_vec.y() * seg_vec.y();
        FT proj = (v.x() * seg_vec.x() + v.y() * seg_vec.y()) / seg_len_sq;

        if (proj > 0 && proj < 1) {
            // Closest point is interior to segment
            Point_2 closest = seg.source() + proj * seg_vec;
            return distance_squared(center, closest) <= radius_sq;
        }

        return false;
    }

    // ========================================================================
    // HELPER FUNCTIONS
    // ========================================================================

    FT distance_squared(const Point_2& p, const Point_2& q) const {
        FT dx = p.x() - q.x();
        FT dy = p.y() - q.y();
        return dx * dx + dy * dy;
    }
};

}  // namespace CGAL

#endif  // CGAL_CIRCLE_SEGMENT_RAY_LINE_EXTENSION_HPP
