#ifndef AABBCOLLISIONDETECTOR_HEADER
#define AABBCOLLISIONDETECTOR_HEADER

#include <CGAL/Minkowski_sum_2/new/aabb/AABB_tree_2.h>
#include <CGAL/Minkowski_sum_2/new/aabb/AABB_traits_2.h>
#include <CGAL/Minkowski_sum_2/new/aabb/AABB_segment_2_primitive.h>

namespace CGAL {

template <class Kernel_, class Container_> class AABB_collision_detector_2 {

public:

    typedef typename Kernel_::Point_2 Point_2;
    typedef typename Kernel_::Vector_2 Vector_2;
    typedef typename CGAL::Polygon_2<Kernel_> Polygon_2;
    typedef typename Polygon_2::Edge_const_iterator Edge_iterator;
    typedef AABB_segment_2_primitive<Kernel_, Edge_iterator, Polygon_2> Tree_segment_2;
    typedef AABB_traits_2<Kernel_, Tree_segment_2> Tree_traits;
    typedef AABB_tree_2<Tree_traits> Tree_2;

public:

    AABB_collision_detector_2(const Polygon_2 &p, const Polygon_2 &q)
        : m_stationary_tree((p.edges_begin()), (p.edges_end())), m_translating_tree((q.edges_begin()), (q.edges_end())), m_p(q), m_q(p) {
    }

    bool check_collision(const Point_2 &t) {
        if (m_stationary_tree.do_intersect_join(m_translating_tree, t, m_p, m_q)) {
            return true;
        }

        Polygon_2 translated_p = transform(typename Kernel_::Aff_transformation_2(CGAL::Translation(), Vector_2(CGAL::ORIGIN, t)), m_p);
        return (translated_p.has_on_bounded_side(*(m_q.vertices_begin())) || m_q.has_on_bounded_side(*(translated_p.vertices_begin())));
    }

private:

    Tree_2 m_stationary_tree;
    Tree_2 m_translating_tree;
    const Polygon_2 &m_p;
    const Polygon_2 &m_q;
};

} // namespace CGAL

#endif
