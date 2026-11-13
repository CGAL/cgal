// ============================================================================
// TASK 4: GEODESIC ARC TRAITS EXTENSION
// Add support for AosLandmarkTraits_2 (event-aware landmark handling)
//
// Geodesic arcs: Curves on sphere/ellipsoid (generalized great circles)
// Landmarks: Poles, discontinuities, cut points on the surface
// ============================================================================

#ifndef CGAL_GEODESIC_ARC_TRAITS_EXTENSION_HPP
#define CGAL_GEODESIC_ARC_TRAITS_EXTENSION_HPP

#include <vector>
#include <utility>

namespace CGAL {

// ============================================================================
// GEODESIC ARC TYPE (on Sphere)
// ============================================================================

template <typename Kernel_>
class Geodesic_arc_3 {
public:
    typedef Kernel_ Kernel;
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Sphere_3 Sphere_3;

private:
    Sphere_3 m_surface;   // Sphere the arc lies on
    Point_3 m_source;     // Starting point
    Point_3 m_target;     // Ending point
    bool m_short_arc;     // true = short arc, false = long arc

public:
    Geodesic_arc_3() : m_short_arc(true) {}

    Geodesic_arc_3(const Sphere_3& srf, const Point_3& s, const Point_3& t,
                   bool short_arc = true)
        : m_surface(srf), m_source(s), m_target(t), m_short_arc(short_arc) {}

    const Sphere_3& surface() const { return m_surface; }
    const Point_3& source() const { return m_source; }
    const Point_3& target() const { return m_target; }
    bool is_short_arc() const { return m_short_arc; }
};

// ============================================================================
// LANDMARK TRAITS FOR GEODESIC ARCS
// ============================================================================
// Handle special points: poles, cut points, discontinuities

template <typename Kernel_>
class Geodesic_arc_landmark_traits_2 {
public:
    typedef Kernel_ Kernel;
    typedef typename Kernel::Point_3 Point_3;
    typedef typename Kernel::Sphere_3 Sphere_3;
    typedef typename Kernel::FT FT;
    typedef Geodesic_arc_3<Kernel> Curve_3;

    // ========================================================================
    // LANDMARK TYPES
    // ========================================================================

    enum Landmark_type {
        POLE,                    // North/South pole on sphere
        CUT_POINT,              // Antipodal point (geodesic discontinuity)
        DISCONTINUITY,          // Other special point
        REGULAR_POINT           // Normal point on arc
    };

    struct Landmark {
        Point_3 point;
        Landmark_type type;
        FT parameter;  // Parametric position along arc [0, 1]
    };

    // ========================================================================
    // IDENTIFY POLES ON SPHERE
    // ========================================================================

    std::vector<Landmark> get_pole_landmarks(const Sphere_3& sphere) const {
        std::vector<Landmark> landmarks;

        // Get sphere center and radius
        const Point_3& center = sphere.center();
        FT radius = sqrt(sphere.squared_radius());

        // North pole (center + radius * Z-axis)
        {
            Landmark north;
            north.point = Point_3(center.x(), center.y(), center.z() + radius);
            north.type = POLE;
            north.parameter = 0.0;  // Parametric position
            landmarks.push_back(north);
        }

        // South pole (center - radius * Z-axis)
        {
            Landmark south;
            south.point = Point_3(center.x(), center.y(), center.z() - radius);
            south.type = POLE;
            south.parameter = 1.0;
            landmarks.push_back(south);
        }

        return landmarks;
    }

    // ========================================================================
    // IDENTIFY ANTIPODAL POINTS (CUT POINTS)
    // ========================================================================
    // Antipodal point creates a *discontinuity* in the geodesic distance

    bool is_antipodal_to(const Point_3& p1, const Point_3& p2,
                        const Sphere_3& sphere) const {
        const Point_3& center = sphere.center();

        // p2 is antipodal to p1 if: center + (center - p1) = p2
        FT tol = 1e-10;

        FT dx = p2.x() - (2 * center.x() - p1.x());
        FT dy = p2.y() - (2 * center.y() - p1.y());
        FT dz = p2.z() - (2 * center.z() - p1.z());

        FT dist_sq = dx * dx + dy * dy + dz * dz;
        return dist_sq < tol * tol;
    }

    // ========================================================================
    // GET ALL LANDMARKS ON GEODESIC ARC
    // ========================================================================

    std::vector<Landmark> get_arc_landmarks(const Curve_3& arc) const {
        std::vector<Landmark> landmarks;

        const Point_3& source = arc.source();
        const Point_3& target = arc.target();
        const Sphere_3& sphere = arc.surface();

        // Poles
        std::vector<Landmark> pole_lms = get_pole_landmarks(sphere);
        for (const auto& p : pole_lms) {
            // Check if pole is on this arc
            if (is_on_arc(arc, p.point)) {
                landmarks.push_back(p);
            }
        }

        // Antipodal point (cut point)
        Point_3 antipodal = get_antipodal_point(target, sphere);
        if (is_on_arc(arc, antipodal)) {
            Landmark cut;
            cut.point = antipodal;
            cut.type = CUT_POINT;
            cut.parameter = 0.5;  // Approximate
            landmarks.push_back(cut);
        }

        return landmarks;
    }

    // ========================================================================
    // LANDMARK-AWARE ARC SPLITTING
    // ========================================================================

    std::vector<Curve_3> split_by_landmarks(const Curve_3& arc) const {
        std::vector<Curve_3> segments;

        std::vector<Landmark> landmarks = get_arc_landmarks(arc);

        // If no internal landmarks, return original arc
        if (landmarks.empty()) {
            segments.push_back(arc);
            return segments;
        }

        // Sort landmarks by parameter
        std::sort(landmarks.begin(), landmarks.end(),
                  [](const Landmark& a, const Landmark& b) {
                      return a.parameter < b.parameter;
                  });

        // Split arc at each landmark
        Point_3 current_source = arc.source();

        for (const auto& lm : landmarks) {
            Curve_3 segment(arc.surface(), current_source, lm.point, true);
            segments.push_back(segment);
            current_source = lm.point;
        }

        // Final segment to target
        {
            Curve_3 final_segment(arc.surface(), current_source, arc.target(), true);
            segments.push_back(final_segment);
        }

        return segments;
    }

    // ========================================================================
    // LANDMARK-AWARE PROPERTIES
    // ========================================================================

    bool arc_passes_through_pole(const Curve_3& arc) const {
        std::vector<Landmark> lms = get_arc_landmarks(arc);
        for (const auto& lm : lms) {
            if (lm.type == POLE) {
                return true;
            }
        }
        return false;
    }

    bool arc_passes_through_cut_point(const Curve_3& arc) const {
        std::vector<Landmark> lms = get_arc_landmarks(arc);
        for (const auto& lm : lms) {
            if (lm.type == CUT_POINT) {
                return true;
            }
        }
        return false;
    }

    // ========================================================================
    // COMPUTE ARC LENGTH (accounting for landmarks)
    // ========================================================================

    FT ComputeArcLength(const Curve_3& arc) const {
        const Point_3& source = arc.source();
        const Point_3& target = arc.target();
        const Sphere_3& sphere = arc.surface();

        FT radius = sqrt(sphere.squared_radius());

        // Compute angle between source and target vectors from center
        const Point_3& center = sphere.center();

        FT dx1 = source.x() - center.x();
        FT dy1 = source.y() - center.y();
        FT dz1 = source.z() - center.z();

        FT dx2 = target.x() - center.x();
        FT dy2 = target.y() - center.y();
        FT dz2 = target.z() - center.z();

        // Dot product
        FT dot = dx1 * dx2 + dy1 * dy2 + dz1 * dz2;
        FT cos_angle = dot / (radius * radius);

        // Clamp to [-1, 1] for numerical stability
        if (cos_angle > 1.0) cos_angle = 1.0;
        if (cos_angle < -1.0) cos_angle = -1.0;

        FT angle = acos(cos_angle);

        // Choose short or long arc
        if (arc.is_short_arc()) {
            return radius * angle;
        } else {
            return radius * (2 * M_PI - angle);
        }
    }

private:
    // ========================================================================
    // HELPER FUNCTIONS
    // ========================================================================

    bool is_on_arc(const Curve_3& arc, const Point_3& p) const {
        // Simplified: check if point is on the great circle
        // In reality, verify: (1) on sphere, (2) between source and target
        const Sphere_3& sphere = arc.surface();

        FT dx = p.x() - sphere.center().x();
        FT dy = p.y() - sphere.center().y();
        FT dz = p.z() - sphere.center().z();

        FT dist_sq = dx * dx + dy * dy + dz * dz;
        FT radius_sq = sphere.squared_radius();

        FT tol = 1e-10;
        return std::abs(dist_sq - radius_sq) < tol;
    }

    Point_3 get_antipodal_point(const Point_3& p, const Sphere_3& sphere) const {
        const Point_3& center = sphere.center();
        return Point_3(2 * center.x() - p.x(),
                       2 * center.y() - p.y(),
                       2 * center.z() - p.z());
    }

    FT acos(FT x) const {
        // Standard acos function
        return std::acos(static_cast<double>(x));
    }

    FT sqrt(FT x) const {
        return std::sqrt(static_cast<double>(x));
    }

    static const FT M_PI;
};

template <typename Kernel_>
const typename Kernel_::FT Geodesic_arc_landmark_traits_2<Kernel_>::M_PI =
    3.14159265358979323846;

}  // namespace CGAL

#endif  // CGAL_GEODESIC_ARC_TRAITS_EXTENSION_HPP
