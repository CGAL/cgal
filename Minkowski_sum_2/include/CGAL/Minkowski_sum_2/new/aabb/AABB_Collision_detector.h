#ifndef AABBCOLLISIONDETECTOR_HEADER
#define AABBCOLLISIONDETECTOR_HEADER

#include "AABB_2d_traits.h"
#include "AABB_segment_2_primitive.h"
#include "AABB_tree_mod.h"

namespace CGAL {
namespace internal {

template <class Kernel_, class Container_> class AABBCollisionDetector {

public:

    typedef typename Kernel_::Point_2 Point;
    typedef typename CGAL::Polygon_2<Kernel_> Polygon_2;
    typedef typename Polygon_2::Traits::Segment_2 Segment_2 ;
    typedef typename Polygon_2::Edge_const_iterator Edge_iterator;
    typedef typename Polygon_2::Edge_const_circulator Edge_circulator;
    typedef AABB_segment_2_primitive<Kernel_, Edge_iterator, Polygon_2> Tree_Segment_2;
    typedef AABB_traits_2<Kernel_, Tree_Segment_2> Tree_Traits;
    typedef AABB_tree<Tree_Traits> AABB_Tree;
    typedef CGAL::Arr_segment_traits_2<Kernel_> Traits_2;

protected:

    Traits_2 m_traits;

public:

    AABBCollisionDetector(Polygon_2 &p, Polygon_2 &q)
        : m_stationary_tree((p.edges_begin()), (p.edges_end())), m_translating_tree((q.edges_begin()), (q.edges_end())), m_p(q), m_q(p) {
    }
    bool checkCollision(const Polygon_2 &p, const Polygon_2 &q) {
        if (m_stationary_tree.do_intersect_join(m_translating_tree, m_translation_point, m_p, m_q)) {
            return true;
        }

        return (p.has_on_bounded_side(*(q.vertices_begin())) || q.has_on_bounded_side(*(p.vertices_begin())));
    }

    void setTranslationPoint(const Point &t) {
        m_translation_point = t;
    }

private:

    AABB_Tree m_stationary_tree;
    AABB_Tree m_translating_tree;
    Point m_translation_point;
    Polygon_2 &m_p;
    Polygon_2 &m_q;
};

} // namespace internal
} // namespace CGAL

#endif
