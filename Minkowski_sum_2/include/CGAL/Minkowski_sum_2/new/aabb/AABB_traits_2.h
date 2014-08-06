#ifndef CGAL_AABB_TRAITS_2_H
#define CGAL_AABB_TRAITS_2_H

namespace CGAL {

double eps(double x) {
    return abs(nextafter(x, DBL_MAX) - x);
}

template<typename GeomTraits, typename AABB_primitive_>
class AABB_traits_2 {
public:
    typedef AABB_traits_2<GeomTraits, AABB_primitive_> AT;
    typedef typename CGAL::Bbox_2 Bounding_box;

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
    typedef typename GeomTraits::Construct_cartesian_const_iterator_2 Construct_cartesian_const_iterator_3;

    AABB_traits_2(const Point &point, const Container &p, const Container &q): m_t_point(point), m_p(p), m_q(q) {
        m_x_interval = Interval_nt<true>(to_interval(point.x()));
        m_y_interval = Interval_nt<true>(to_interval(point.y()));
        m_px = to_double(point.x());
        m_py = to_double(point.y());
    };

    AABB_traits_2(): m_p(Container()), m_q(Container()) {
    };

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

            if (bbox.xmax()-bbox.xmin() >= bbox.ymax()-bbox.ymin()) {
                std::nth_element(first, middle, beyond, less_x); // sort along x
            } else {
                std::nth_element(first, middle, beyond, less_y); // sort along y
            }
        }
    };

    Sort_primitives sort_primitives_object() const {
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
            typename AT::Bounding_box bbox = first->datum().bbox();

            for (++first; first != beyond; ++first) {
                bbox = bbox + first->datum().bbox();
            }

            return bbox;
        }
    };

    Compute_bbox compute_bbox_object() const {
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
            double x_epsilon = max(max(eps(m_traits->m_px),eps(bbox.xmin())),eps(bbox.xmax()))*2;
            // Get y max error.
            double y_epsilon = max(max(eps(m_traits->m_py),eps(bbox.ymin())),eps(bbox.ymax()))*2;
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

            return do_overlap(q, t_box);
        }

        bool operator()(const Primitive &q, const Bounding_box &bbox) const {
            /* Code for faster bbox, needs to be tested
            // Get x max error.
            double x_epsilon = max(max(eps(m_traits->m_px),eps(bbox.xmin())),eps(bbox.xmax()))*2;
            // Get y max error.
            double y_epsilon = max(max(eps(m_traits->m_py),eps(bbox.ymin())),eps(bbox.ymax()))*2;
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

            return do_overlap(q.datum().bbox(), t_box);
        }

        bool operator()(const Bounding_box &q, const Primitive &pr) const {

            typename Primitive::Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(Translation(), Vector_2(ORIGIN, m_traits->get_translation_point())));
            return do_overlap(q, tr_pr.bbox());
        }

        bool operator()(const Primitive &q, const Primitive &pr) const {

            typename Primitive::Datum tr_pr = pr.datum().transform(typename GeomTraits::Aff_transformation_2(Translation(), Vector_2(ORIGIN, m_traits->get_translation_point())));

            if (!do_overlap(q.datum().bbox(), tr_pr.bbox())) {
                return false;
            }

            return GeomTraits().intersect_2_object()(q.datum(), tr_pr);
        }
    };

    Do_intersect do_intersect_object() {
        return Do_intersect(this);
    }

private:

    Point m_t_point;
    double m_px, m_py;
    Interval_nt<true> m_x_interval;
    Interval_nt<true> m_y_interval;
    const Container &m_p;
    const Container &m_q;

    /// Comparison functions
    static bool less_x(const Primitive &pr1, const Primitive &pr2) {
        return pr1.reference_point().x() < pr2.reference_point().x();
    }

    static bool less_y(const Primitive &pr1, const Primitive &pr2) {
        return pr1.reference_point().y() < pr2.reference_point().y();
    }
};

} // namespace CGAL

#endif
