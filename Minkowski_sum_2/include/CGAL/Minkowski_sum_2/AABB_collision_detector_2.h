#ifndef CGAL_AABB_COLLISION_DETECTOR_2_H
#define CGAL_AABB_COLLISION_DETECTOR_2_H

#include <CGAL/Minkowski_sum_2/AABB_tree_with_join.h>
#include <CGAL/Minkowski_sum_2/AABB_traits_2.h>
#include <CGAL/Minkowski_sum_2/AABB_segment_2_primitive.h>

namespace CGAL {

// Tests whether two polygons P and Q overlap for different translations of Q.
template <class Kernel_, class Container_>
class AABB_collision_detector_2
{

public:

  typedef typename Kernel_::Point_2 Point_2;
  typedef typename Kernel_::Vector_2 Vector_2;
  typedef typename CGAL::Polygon_2<Kernel_> Polygon_2;
  typedef typename Polygon_2::Edge_const_iterator Edge_iterator;
  typedef AABB_segment_2_primitive<Kernel_, Edge_iterator, Polygon_2>
  Tree_segment_2;
  typedef AABB_traits_2<Kernel_, Tree_segment_2> Tree_traits;
  typedef AABB_tree_with_join<Tree_traits> Tree_2;

public:

  AABB_collision_detector_2(const Polygon_2 &p, const Polygon_2 &q)
    : m_stationary_tree((p.edges_begin()), (p.edges_end())),
      m_translating_tree((q.edges_begin()), (q.edges_end())), m_p(q), m_q(p)
  {
  }

  // Returns true iff the polygons' boundaries intersect or one polygon is
  // completely inside of the other one. Q is translated by t.
  bool check_collision(const Point_2 &t)
  {
    if (m_stationary_tree.do_intersect(m_translating_tree, t))
    {
      return true;
    }

    // If t_q is inside of P, or t_p is inside of Q, one polygon is completely
    // inside of the other.
    Point_2 t_q = *m_q.vertices_begin() + Vector_2(ORIGIN, t);
    Point_2 t_p = *m_p.vertices_begin() - Vector_2(ORIGIN, t);

    // Use bounded_side_2() instead of on_bounded_side() because the latter
    // checks vor simplicity every time.
    return bounded_side_2(m_p.vertices_begin(), m_p.vertices_end(), t_q, m_p.traits_member()) == ON_BOUNDED_SIDE
        || bounded_side_2(m_q.vertices_begin(), m_q.vertices_end(), t_p, m_q.traits_member()) == ON_BOUNDED_SIDE;
  }

private:

  Tree_2 m_stationary_tree;
  Tree_2 m_translating_tree;
  const Polygon_2 &m_p;
  const Polygon_2 &m_q;
};

} // namespace CGAL

#endif
