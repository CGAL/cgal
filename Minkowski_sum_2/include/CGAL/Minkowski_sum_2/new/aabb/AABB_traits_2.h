#ifndef CGAL_AABB_2D_TRAITS_H_
#define CGAL_AABB_2D_TRAITS_H_

namespace CGAL {
namespace internal {

template <class R_> class Point_2;
template <class R_> class Segment_2;

double eps(double x) {
    return CGAL::abs(CGAL::nextafter(x, DBL_MAX) - x);
}

template<typename GeomTraits, typename AABB_primitive_>
class AABB_traits_2 {
public:
    typedef AABB_traits_2<GeomTraits, AABB_primitive_> AT;
    typedef typename CGAL::Bbox_2 Bounding_box;
    typedef typename CGAL::Object Object;

    typedef AABB_primitive_ Primitive;
    typedef typename Primitive::Id Id;
    typedef typename Primitive::Datum Datum;
    typedef typename Primitive::Container Container;

    typedef typename GeomTraits::Point_2 Point;
    typedef typename GeomTraits::Vector_2 Vector_2;

    typedef typename std::pair<Object, typename Primitive::Id> Object_and_primitive_id;
    typedef typename std::pair<Point, typename Primitive::Id> Point_and_primitive_id;

    // types for search tree
    typedef typename GeomTraits::FT FT;
    typedef typename GeomTraits::Point_2 Point_3;
    typedef typename GeomTraits::Circle_2 Sphere_3;
    typedef typename GeomTraits::Iso_rectangle_2 Iso_cuboid_3;
    typedef typename GeomTraits::Construct_center_2 Construct_center_3;
    typedef typename GeomTraits::Construct_iso_rectangle_2 Construct_iso_cuboid_3;
    typedef typename GeomTraits::Construct_min_vertex_2 Construct_min_vertex_3;
    typedef typename GeomTraits::Construct_max_vertex_2 Construct_max_vertex_3;
    typedef typename GeomTraits::Compute_squared_radius_2 Compute_squared_radius_3;
    typedef typename GeomTraits::Compute_squared_distance_2 Compute_squared_distance_3;
    typedef typename GeomTraits::Cartesian_const_iterator_2 Cartesian_const_iterator_3;
    typedef typename GeomTraits::Construct_cartesian_const_iterator_2
    Construct_cartesian_const_iterator_3;

    AABB_traits_2(const Point &point, const Container &p, const Container &q): m_t_point(point), m_p(p), m_q(q) {
        m_x_interval = Interval_nt<true>(CGAL::to_interval(point.x()));
        m_y_interval = Interval_nt<true>(CGAL::to_interval(point.y()));
        m_px = CGAL::to_double(point.x());
        m_py = CGAL::to_double(point.y());
    };

    AABB_traits_2(): m_p(Container()), m_q(Container()) {
    };

    ~AABB_traits_2() { };

    Interval_nt<true> get_int_x() const {
        return m_x_interval;
    }
    Interval_nt<true> get_int_y() const {
        return m_y_interval;
    }

    Point get_translation_point() const {
        return m_t_point;
    }
    const Container &get_p() const {
        return m_p;
    }
    const Container &get_q() const {
        return m_q;
    }

    /**
     * @brief Sorts [first,beyond[
     * @param first iterator on first element
     * @param beyond iterator on beyond element
     * @param bbox the bounding box of [first,beyond[
     *
     * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longuest
     * axis, using the comparison function <dim>_less_than (dim in {x,y,z})
     */
    class Sort_primitives {

    public:

        template<typename PrimitiveIterator>
        void operator()(PrimitiveIterator first,
                        PrimitiveIterator beyond,
                        const typename AT::Bounding_box &bbox) const {
            PrimitiveIterator middle = first + (beyond - first) / 2;

            switch (longest_axis(bbox)) {
            case AT::CGAL_AXIS_X: // sort along x
                std::nth_element(first, middle, beyond, less_x);
                break;

            case AT::CGAL_AXIS_Y: // sort along y
                std::nth_element(first, middle, beyond, less_y);
                break;

            case AT::CGAL_AXIS_Z: // sort along z
                CGAL_error();
                break;

            default:
                CGAL_error();
            }
        }
    };

    Sort_primitives sort_primitives_object() {
        return Sort_primitives();
    }

    /**
     * Computes the bounding box of a set of primitives
     * @param first an iterator on the first primitive
     * @param beyond an iterator on the past-the-end primitive
     * @return the bounding box of the primitives of the iterator range
     */
    class Compute_bbox {

    public:

        template<typename ConstPrimitiveIterator>
        typename AT::Bounding_box operator()(ConstPrimitiveIterator first,
                                             ConstPrimitiveIterator beyond) const {
            typename AT::Bounding_box bbox = compute_bbox(*first);

            for (++first; first != beyond; ++first) {
                bbox = bbox + compute_bbox(*first);
            }

            return bbox;
        }
    };

    Compute_bbox compute_bbox_object() {
        return Compute_bbox();
    }

    class Do_intersect {

    private:

        AABB_traits_2 *m_traits;
        typedef typename Primitive::Datum Datum;

    public:

        Do_intersect(AABB_traits_2 *_traits): m_traits(_traits) {}

        bool operator()(const Bounding_box &q, const Bounding_box &bbox) const {

            /* Code for faster bbox, needs to be tested
            // Get x max error.
            double x_epsilon = CGAL::max(CGAL::max(eps(m_traits->m_px),eps(bbox.xmin())),eps(bbox.xmax()))*2;
            // Get y max error.
            double y_epsilon = CGAL::max(CGAL::max(eps(m_traits->m_py),eps(bbox.ymin())),eps(bbox.ymax()))*2;
            double t_left = (m_traits->m_px + bbox.xmin())-x_epsilon;
            double t_right = (m_traits->m_px + bbox.xmax())+x_epsilon;
            double t_bottom = (m_traits->m_py + bbox.ymin())-y_epsilon;
            double t_top = (m_traits->m_py + bbox.ymax())+y_epsilon;
            Bounding_box t_box(t_left,t_bottom,t_right,t_top);
            */

            double t_left = (m_traits->get_int_x() + bbox.xmin()).inf();
            double t_right = (m_traits->get_int_x() + bbox.xmax()).sup();
            double t_bottom = (m_traits->get_int_y() + bbox.ymin()).inf();
            double t_top = (m_traits->get_int_y() + bbox.ymax()).sup();
            Bounding_box t_box(t_left, t_bottom, t_right, t_top);

            return CGAL::do_overlap(q, t_box);
        }

        bool operator()(const Primitive &q, const Bounding_box &bbox) const {
            /* Code for faster bbox, needs to be tested
            // Get x max error.
            double x_epsilon = CGAL::max(CGAL::max(eps(m_traits->m_px),eps(bbox.xmin())),eps(bbox.xmax()))*2;
            // Get y max error.
            double y_epsilon = CGAL::max(CGAL::max(eps(m_traits->m_py),eps(bbox.ymin())),eps(bbox.ymax()))*2;
            double t_left = (m_traits->m_px + bbox.xmin())-x_epsilon;
            double t_right = (m_traits->m_px + bbox.xmax())+x_epsilon;
            double t_bottom = (m_traits->m_py + bbox.ymin())-y_epsilon;
            double t_top = (m_traits->m_py + bbox.ymax())+y_epsilon;
            Bounding_box t_box(t_left,t_bottom,t_right,t_top);
            */

            double t_left = (m_traits->get_int_x() + bbox.xmin()).inf();
            double t_right = (m_traits->get_int_x() + bbox.xmax()).sup();
            double t_bottom = (m_traits->get_int_y() + bbox.ymin()).inf();
            double t_top = (m_traits->get_int_y() + bbox.ymax()).sup();
            Bounding_box t_box(t_left, t_bottom, t_right, t_top);

            return CGAL::do_overlap(q.datum().bbox(), t_box);
        }

        bool operator()(const Bounding_box &q, const Primitive &pr) const {

            typename Primitive::Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, m_traits->get_translation_point())));
            return CGAL::do_overlap(q, tr_pr.bbox());
        }

        bool operator()(const Primitive &q, const Primitive &pr) const {

            typename Primitive::Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, m_traits->get_translation_point())));

            if (!CGAL::do_overlap(q.datum().bbox(), tr_pr.bbox())) {
                return false;
            }

            CGAL::Object intersection_object = GeomTraits().intersect_2_object()(q.datum(), tr_pr);

            if (const CGAL::Point_2<GeomTraits> *ipoint = CGAL::object_cast<CGAL::Point_2<GeomTraits> >(&intersection_object)) {
                // handle weak intersections
                bool has_weak_intersection = false;
                bool p_intersect = false;
                bool p_intersect_start = false;
                bool q_intersect = false;
                bool q_intersect_start = false;

                if (*ipoint == tr_pr.source()) {
                    has_weak_intersection = true;
                    p_intersect = true;
                    p_intersect_start = true;
                } else {
                    if (*ipoint == tr_pr.target()) {
                        has_weak_intersection = true;
                        p_intersect = true;
                    }
                }

                if (*ipoint == q.datum().source()) {
                    has_weak_intersection = true;
                    q_intersect = true;
                    q_intersect_start = true;
                } else {
                    if (*ipoint == q.datum().target()) {
                        q_intersect = true;
                        has_weak_intersection = true;
                    }
                }

                if (has_weak_intersection) {

                    bool val = handle_weak_intersections(p_intersect, q_intersect, p_intersect_start, q_intersect_start, pr, q, tr_pr);

                    if (val == false) {
                        int k = 4;
                        k = k + 4;
                        k++;
                    }

                    return val;
                } else {
                    return true;
                }
            }

            if (const CGAL::Segment_2<GeomTraits> *iseg = CGAL::object_cast<CGAL::Segment_2<GeomTraits> >(&intersection_object)) { // we have overlapping segments
                CGAL::Comparison_result c1 = CGAL::compare_xy(tr_pr.source(), tr_pr.target());
                CGAL::Comparison_result c2 = CGAL::compare_xy(q.datum().source(), q.datum().target());

                bool same_dir = (c1 == c2);
                return same_dir;
            } else {
                return false; // no intersection
            }
        }

    private:

        bool handle_weak_intersections(bool p_intersect, bool q_intersect, bool p_intersect_start, bool q_intersect_start, const Primitive &p, const Primitive &q, const Datum &tr_pr_datum) const {
            Id itr_p = p.id();
            Id itr_q = q.id();
            Id p_other = get_other_segment(p_intersect_start, itr_p, m_traits->get_p());
            Id q_other = get_other_segment(q_intersect_start, itr_q, m_traits->get_q());
            Datum p_other_translated = (*p_other).transform(typename GeomTraits::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, m_traits->get_translation_point())));

            if (p_intersect && !q_intersect) {
                if (p_intersect_start) {
                    return handle_weak_intersection(p_other_translated, tr_pr_datum, *itr_q);
                } else {
                    return handle_weak_intersection(tr_pr_datum, p_other_translated, *itr_q);
                }
            } else {
                if (!p_intersect && q_intersect) {
                    if (q_intersect_start) {
                        return handle_weak_intersection(*q_other, *itr_q, tr_pr_datum);
                    } else {
                        return handle_weak_intersection(*itr_q, *q_other, tr_pr_datum);
                    }
                } else {
                    Datum first_p, second_p;
                    Datum first_q, second_q;

                    if (p_intersect_start) {
                        first_p = p_other_translated;
                        second_p = tr_pr_datum;
                    } else {
                        first_p = tr_pr_datum;
                        second_p = p_other_translated;
                    }

                    if (q_intersect_start) {
                        first_q = *q_other;
                        second_q = *itr_q;
                    } else {
                        first_q = *itr_q;
                        second_q = *q_other;
                    }

                    return is_overlapping(first_p, second_p, first_q, second_q);
                }
            }
        }

        bool handle_weak_intersection(const Datum &incoming, const Datum &outgoing, const Datum &other_segment) const {
            // There is an overlap in polygon regions if the outgoing of p is ccw-between outgoing q and -incoming q or vice versa.
            //return (other_segment.direction()).counterclockwise_in_between(outgoing.direction(),incoming.opposite().direction());
            return (other_segment.direction()).counterclockwise_in_between(outgoing.direction(), incoming.opposite().direction()) ||
                   outgoing.direction().counterclockwise_in_between(other_segment.direction(), other_segment.opposite().direction());
        }

        bool is_overlapping(const Datum &incoming_p, const Datum &outgoing_p, const Datum &incoming_q, const Datum &outgoing_q) const {
            // There is an overlap in polygon regions if the outgoing of p is ccw-between outgoing q and -incoming q or vice versa.
            return ((outgoing_q.direction()).counterclockwise_in_between(outgoing_p.direction(), incoming_p.opposite().direction()) ||
                    (outgoing_p.direction()).counterclockwise_in_between(outgoing_q.direction(), incoming_q.opposite().direction()));
        }

        Id get_other_segment(bool start, const Id &itr_p, const Container &cont) const {
            Id p_other;

            if (start) {
                p_other = cont.edges_begin();

                if (p_other == itr_p) {
                    p_other = cont.edges_end();
                    --p_other;
                } else {
                    while (p_other != itr_p) {
                        ++p_other;
                    }

                    --p_other;
                }
            } else {
                p_other = cont.edges_end();
                --p_other;

                if (p_other == itr_p) {
                    p_other = cont.edges_begin();
                } else {
                    while (p_other != itr_p) {
                        --p_other;
                    }

                    ++p_other;
                }
            }

            return p_other;
        }
    };

    Do_intersect do_intersect_object() {
        return Do_intersect(this);
    }

    class Intersection {

    public:

        template<typename Query>
        boost::optional<typename AT::Object_and_primitive_id>
        operator()(const Query &query, const typename AT::Primitive &primitive) const {
            typedef boost::optional<Object_and_primitive_id> Intersection;

            CGAL::Object object = GeomTraits().intersect_2_object()(primitive.datum(), query);

            if (object.empty()) {
                return Intersection();
            } else {
                return Intersection(Object_and_primitive_id(object, primitive.id()));
            }
        }
    };

    Intersection intersection_object() {
        return Do_intersect(this);
    }

    // This should go down to the GeomTraits, i.e. the kernel
    class Closest_point {

        typedef typename AT::Point Point;
        typedef typename AT::Primitive Primitive;

    public:

        Point operator()(const Point &p, const Primitive &pr, const Point &bound) const {
            // seems to be unused:
            //return CGAL::nearest_point_2(p, pr.datum(), bound);
            return p;
        }
    };

    // This should go down to the GeomTraits, i.e. the kernel
    // and the internal implementation should change its name from
    // do_intersect to something like does_contain (this is what we compute,
    // this is not the same do_intersect as the spherical kernel)
    class Compare_distance {

        typedef typename AT::Point Point;
        typedef typename AT::Primitive Primitive;

    public:

        template <class Solid>
        CGAL::Comparison_result operator()(const Point &p, const Solid &pr, const Point &bound) const {
            return GeomTraits().do_intersect_2_object()
                   (GeomTraits().construct_sphere_2_object()
                    (p, GeomTraits().compute_squared_distance_2_object()(p, bound)), pr) ?
                   CGAL::SMALLER : CGAL::LARGER;
        }
    };

    Closest_point closest_point_object() {
        return Closest_point();
    }

    Compare_distance compare_distance_object() {
        return Compare_distance();
    }

private:

    Point m_t_point;
    double m_px, m_py;
    Interval_nt<true> m_x_interval;
    Interval_nt<true> m_y_interval;
    const Container &m_p;
    const Container &m_q;

    /**
     * @brief Computes bounding box of one primitive
     * @param pr the primitive
     * @return the bounding box of the primitive \c pr
     */
    static Bounding_box compute_bbox(const Primitive &pr) {
        return pr.datum().bbox();
    }

    typedef enum {
        CGAL_AXIS_X = 0,
        CGAL_AXIS_Y = 1,
        CGAL_AXIS_Z = 2
    } Axis;

    static Axis longest_axis(const Bounding_box &bbox);

    /// Comparison functions
    static bool less_x(const Primitive &pr1, const Primitive &pr2) {
        return pr1.reference_point().x() < pr2.reference_point().x();
    }

    static bool less_y(const Primitive &pr1, const Primitive &pr2) {
        return pr1.reference_point().y() < pr2.reference_point().y();
    }
};

template<typename GT, typename P>
typename AABB_traits_2<GT, P>::Axis
AABB_traits_2<GT, P>::longest_axis(const Bounding_box &bbox) {
    const double dx = bbox.xmax() - bbox.xmin();
    const double dy = bbox.ymax() - bbox.ymin();

    if (dx >= dy) {
        return CGAL_AXIS_X;
    } else {
        return CGAL_AXIS_Y;
    }
}

} // namespace internal
} // namespace CGAL

#endif
