// compute self-intersection of a CGAL triangle polyhedron mesh
// original code from Lutz Kettner

#ifndef CGAL_SELF_INTERSECTION_POLYHEDRON_3
#define CGAL_SELF_INTERSECTION_POLYHEDRON_3

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#include <vector>
#include <exception>

#include <boost/function_output_iterator.hpp>

namespace CGAL {
namespace internal {
template <class Polyhedron, class Kernel, class OutputIterator>
struct Intersect_facets
{
// typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;
  typedef typename Polyhedron::Facet_const_handle       Facet_const_handle;
  typedef typename Polyhedron::Halfedge_const_handle    Halfedge_const_handle;
  typedef typename CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Facet_const_handle> Box;
// members
  mutable OutputIterator m_iterator;
  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_facets(OutputIterator it, const Kernel& kernel)
    : m_iterator(it),
    segment_functor(kernel.construct_segment_3_object()),
    triangle_functor(kernel.construct_triangle_3_object()),
    do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b,
    const Box* c) const
  {
    Halfedge_const_handle h = b->handle()->halfedge();

    // check for shared egde --> no intersection
    if(h->opposite()->facet() == c->handle() ||
      h->next()->opposite()->facet() == c->handle() ||
      h->next()->next()->opposite()->facet() == c->handle())
      return;

    // check for shared vertex --> maybe intersection, maybe not
    Halfedge_const_handle g = c->handle()->halfedge();
    Halfedge_const_handle v;

    if(h->vertex() == g->vertex())
      v = g;
    if(h->vertex() == g->next()->vertex())
      v = g->next();
    if(h->vertex() == g->next()->next()->vertex())
      v = g->next()->next();

    if(v == Halfedge_const_handle())
    {
      h = h->next();
      if(h->vertex() == g->vertex())
  v = g;
      if(h->vertex() == g->next()->vertex())
  v = g->next();
      if(h->vertex() == g->next()->next()->vertex())
  v = g->next()->next();
      if(v == Halfedge_const_handle())
      {
  h = h->next();
  if(h->vertex() == g->vertex())
    v = g;
  if(h->vertex() == g->next()->vertex())
    v = g->next();
  if(h->vertex() == g->next()->next()->vertex())
    v = g->next()->next();
      }
    }

    if(v != Halfedge_const_handle())
    {
      // found shared vertex: 
      CGAL_assertion(h->vertex() == v->vertex());
      // geometric check if the opposite segments intersect the triangles
      Triangle t1 = triangle_functor( h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
      Triangle t2 = triangle_functor( v->vertex()->point(), v->next()->vertex()->point(), v->next()->next()->vertex()->point());

      Segment s1 = segment_functor( h->next()->vertex()->point(), h->next()->next()->vertex()->point());
      Segment s2 = segment_functor( v->next()->vertex()->point(), v->next()->next()->vertex()->point());

      if(do_intersect_3_functor(t1,s2))
      {
        *m_iterator++ = std::make_pair(b->handle(), c->handle());
      }
      else if(do_intersect_3_functor(t2,s1))
      {
        *m_iterator++ = std::make_pair(b->handle(), c->handle());
      }
      return;
    }

    // check for geometric intersection
    Triangle t1 = triangle_functor( h->vertex()->point(), h->next()->vertex()->point(), h->next()->next()->vertex()->point());
    Triangle t2 = triangle_functor( g->vertex()->point(), g->next()->vertex()->point(), g->next()->next()->vertex()->point());
    if(do_intersect_3_functor(t1, t2))
    {
      *m_iterator++ = std::make_pair(b->handle(), c->handle());
    }
  } // end operator ()
}; // end struct Intersect_facets

struct Throw_at_output {
  class Throw_at_output_exception: public std::exception
  { };

  template<class T>
  void operator()(const T& /* t */) const {
    throw Throw_at_output_exception();
  }
};

}// namespace internal

////////////////////////////////////////////////////////////////////////////////////
/*
/// Geometric traits concept for the functions `self_intersect`
concept SelfIntersectionTraits{
  /// @name Geometric Types
  /// @{
  /// 3D point type
  typedef unspecified_type Point_3;
  /// 3D triangle type
  typedef unspecified_type Triangle_3;
  /// 3D segment type
  typedef unspecified_type Segment_3;
  /// @}

  /// @name Functors
  /// @{
  /// Functor constructing triangles. It provides `Triangle_3 operator() const(const Point_3&, const Point_3&, const Point_3&)
  typedef unspecified_type Construct_triangle_3;
  /// Functor constructing segments. It provides `Segment_3 operator() const(const Point_3&, const Point_3&)
  typedef unspecified_type Construct_segment_3;
  /// Functor testing intersections between triangles and segment. It provides `bool operator() const (const Triangle_3&, const Segment_3&)` and `bool operator() const (const Triangle_3&, const Triangle_3&)`
  typedef unspecified_type Do_intersect_3;
  /// @}

  /// @name Functions
  /// @{
  Construct_triangle_3 construct_triangle_3_object() const;
  Construct_segment_3 construct_segment_3_object() const;
  Do_intersect_3 do_intersect_3_object() const;
  /// @}
};
*/
////////////////////////////////////////////////////////////////////////////////////

/** 
 * Detects and reports self-intersections of a triangulated polyhedral surface
 * @pre @a p.is_pure_triangle()
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam Polyhedron a \cgal polyhedron
 * @tparam OutputIterator Output iterator accepting objects of type `std::pair<Polyhedron::Facet_const_handle, Polyhedron::Facet_const_handle>`
 *
 * @param polyhedron polyhedron to be checked
 * @param out all pairs of non-adjacent facets intersecting are put in it
 *
 * @return out
 *
 * \TODO Doc: move SelfIntersectionTraits concept to appropriate location.
 */
template <class GeomTraits, class Polyhedron, class OutputIterator>
OutputIterator self_intersect(const Polyhedron& polyhedron, OutputIterator out, const GeomTraits& geom_traits = GeomTraits())
{
  CGAL_assertion(polyhedron.is_pure_triangle());

  typedef typename Polyhedron::Facet_const_iterator     Facet_const_iterator;
  typedef typename Polyhedron::Facet_const_handle       Facet_const_handle;
  typedef typename CGAL::Box_intersection_d::Box_with_handle_d<double, 3, Facet_const_handle> Box;

  // make one box per facet
  std::vector<Box> boxes;
  boxes.reserve(polyhedron.size_of_facets());

  Facet_const_iterator f;
  for(f = polyhedron.facets_begin();
    f != polyhedron.facets_end();
    f++)
    boxes.push_back(Box( f->halfedge()->vertex()->point().bbox() +
    f->halfedge()->next()->vertex()->point().bbox() +
    f->halfedge()->next()->next()->vertex()->point().bbox(),
    f));

  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(polyhedron.size_of_facets());
  typename std::vector<Box>::iterator b;
  for(b = boxes.begin();
    b != boxes.end();
    b++)
    box_ptr.push_back(&*b);

  // compute self-intersections filtered out by boxes
  internal::Intersect_facets<Polyhedron,GeomTraits,OutputIterator> intersect_facets(out, geom_traits);
  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}

/** 
 * Checks if a polyhedron is self-intersecting
 * @pre @a p.is_pure_triangle()
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam Polyhedron a %CGAL polyhedron
 *
 * @param polyhedron polyhedron to be tested
 *
 * @return true if `polyhedron` is self-intersecting
 */
template <class GeomTraits, class Polyhedron>
bool self_intersect(const Polyhedron& polyhedron, const GeomTraits& geom_traits = GeomTraits())
{
  try 
  {
    typedef boost::function_output_iterator<internal::Throw_at_output> OutputIterator;
    self_intersect<GeomTraits>(polyhedron, OutputIterator(), geom_traits); 
  }
  catch( internal::Throw_at_output::Throw_at_output_exception& ) 
  { return true; }

  return false;
}

}// namespace CGAL

#endif // CGAL_SELF_INTERSECTION_POLYHEDRON_3
