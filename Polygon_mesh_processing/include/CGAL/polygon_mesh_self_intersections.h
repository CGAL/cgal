// Copyright (c) 2008 INRIA Sophia-Antipolis (France).
// Copyright (c) 2008-2013 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Pierre Alliez, Laurent Rineau, Ilker O. Yaz

// compute self-intersection of a CGAL triangle polyhedron mesh
// original code from Lutz Kettner

#ifndef CGAL_POLYGON_MESH_SELF_INTERSECTIONS
#define CGAL_POLYGON_MESH_SELF_INTERSECTIONS

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#include <vector>
#include <exception>

#include <boost/function_output_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>

namespace CGAL {
namespace internal {
template <class TM,//TriangleMesh
          class Kernel,
          class Box, class OutputIterator>
struct Intersect_facets
{
  // wrapper to check whether anything is inserted to output iterator
  struct Output_iterator_with_bool 
  {
    Output_iterator_with_bool(OutputIterator* out, bool* intersected)
      : m_iterator(out), m_intersected(intersected) { }

    template<class T>
    void operator()(const T& t) {
      *m_intersected = true;
      *(*m_iterator)++ = t; 
    }

    OutputIterator* m_iterator;
    bool* m_intersected;
  };
// typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

// members
  const TM& m_tmesh;
  const Ppmap m_point;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  
  Intersect_facets(const TM& tmesh, OutputIterator it, const Kernel& kernel)
    : 
    m_tmesh(tmesh),
    m_point(get(vertex_point, m_tmesh)),
    m_iterator(it),
    m_intersected(false),
    m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
    segment_functor(kernel.construct_segment_3_object()),
    triangle_functor(kernel.construct_triangle_3_object()),
    do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b,
    const Box* c) const
  {
    halfedge_descriptor h  = halfedge(b->info(),m_tmesh);

    // check for shared egde --> no intersection
    if(face(opposite(h,m_tmesh),m_tmesh) == c->info() ||
       face(opposite(next(h,m_tmesh),m_tmesh),m_tmesh) == c->info() ||
       face(opposite(next(next(h,m_tmesh),m_tmesh),m_tmesh),m_tmesh) == c->info())
      return;

    // check for shared vertex --> maybe intersection, maybe not
    halfedge_descriptor g = halfedge(c->info(),m_tmesh);
    halfedge_descriptor v;

    if(target(h,m_tmesh) == target(g,m_tmesh))
      v = g;
    if(target(h,m_tmesh) == target(next(g,m_tmesh),m_tmesh))
      v = next(g,m_tmesh);
    if(target(h,m_tmesh) == target(next(next(g,m_tmesh),m_tmesh),m_tmesh))
      v = next(next(g,m_tmesh),m_tmesh);

    if(v == halfedge_descriptor()){
      h = next(h,m_tmesh);
      if(target(h,m_tmesh) == target(g,m_tmesh))
        v = g;
      if(target(h,m_tmesh) == target(next(g,m_tmesh),m_tmesh))
        v = next(g,m_tmesh);
      if(target(h,m_tmesh) == target(next(next(g,m_tmesh),m_tmesh),m_tmesh))
        v = next(next(g,m_tmesh),m_tmesh);
      if(v == halfedge_descriptor()){
        h = next(h,m_tmesh);
        if(target(h,m_tmesh) == target(g,m_tmesh))
          v = g;
        if(target(h,m_tmesh) == target(next(g,m_tmesh),m_tmesh))
          v = next(g,m_tmesh);
        if(target(h,m_tmesh) == target(next(next(g,m_tmesh),m_tmesh),m_tmesh))
          v = next(next(g,m_tmesh),m_tmesh);
      }
    }

    if(v != halfedge_descriptor()){
      // found shared vertex: 
      CGAL_assertion(target(h,m_tmesh) == target(v,m_tmesh));
      // geometric check if the opposite segments intersect the triangles
      Triangle t1 = triangle_functor( m_point[target(h,m_tmesh)], m_point[target(next(h,m_tmesh),m_tmesh)], m_point[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
      Triangle t2 = triangle_functor( m_point[target(v,m_tmesh)], m_point[target(next(v,m_tmesh),m_tmesh)], m_point[target(next(next(v,m_tmesh),m_tmesh),m_tmesh)]);
      
      Segment s1 = segment_functor( m_point[target(next(h,m_tmesh),m_tmesh)], m_point[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
      Segment s2 = segment_functor( m_point[target(next(v,m_tmesh),m_tmesh)], m_point[target(next(next(v,m_tmesh),m_tmesh),m_tmesh)]);
      
      if(do_intersect_3_functor(t1,s2)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      } else if(do_intersect_3_functor(t2,s1)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      }
      return;
    }
    
    // check for geometric intersection
    Triangle t1 = triangle_functor( m_point[target(h,m_tmesh)], m_point[target(next(h,m_tmesh),m_tmesh)], m_point[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
    Triangle t2 = triangle_functor( m_point[target(g,m_tmesh)], m_point[target(next(g,m_tmesh),m_tmesh)], m_point[target(next(next(g,m_tmesh),m_tmesh),m_tmesh)]);
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
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

namespace Polygon_mesh_processing {

/** 
 * \ingroup PkgPolygonMeshProcessing
 * Detects and reports self-intersections of a triangulated polyhedral surface.
 * Depends on \ref PkgBoxIntersectionDSummary
 * @pre @a CGAL::is_pure_triangle(tmesh)
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam TriangleMesh a model of `FaceListGraph` (possibly a \cgal polyhedron)
 * @tparam OutputIterator Output iterator accepting objects of type 
 *   `std::pair<TriangleMesh::face_descriptor, TriangleMesh::face_descriptor>`
 *   if @a polygon mesh is passed by const reference.
 *
 * @param tmesh triangle mesh to be checked, might be passed by const reference or reference
 * @param out all pairs of non-adjacent facets intersecting are put in it
 * @param geom_traits traits class providing intersection test primitives
 *
 * @return `out`. Note the OutputIterator can be empty.
 */
template <class GeomTraits, class TriangleMesh, class OutputIterator>
OutputIterator
self_intersections(const TriangleMesh& tmesh,
               OutputIterator out,
               const GeomTraits& geom_traits = GeomTraits())
{
  CGAL_precondition(CGAL::is_pure_triangle(tmesh));

  typedef TriangleMesh TM;

  typedef typename boost::graph_traits<TM>::face_iterator Facet_it;

  typedef typename boost::graph_traits<TM>::face_descriptor Facet_hdl;

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, Facet_hdl> Box;

  typedef typename boost::property_map<TM, CGAL::vertex_point_t>::const_type Ppmap;

  Ppmap m_point = get(CGAL::vertex_point, tmesh);

  // make one box per facet
  std::vector<Box> boxes;
  boxes.reserve(num_faces(tmesh));

  Facet_it fi,e;

  for(boost::tie(fi,e)= faces(tmesh);
    fi != e;
      ++fi){
    Facet_hdl f = *fi;
    boxes.push_back(Box( m_point[target(halfedge(f,tmesh),tmesh)].bbox() +
                         m_point[target(next(halfedge(f,tmesh),tmesh),tmesh)].bbox() +
                         m_point[target(next(next(halfedge(f,tmesh),tmesh),tmesh),tmesh)].bbox(),
    f));
  }
  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(num_faces(tmesh));
  typename std::vector<Box>::iterator b;
  for(b = boxes.begin();
    b != boxes.end();
    b++)
    box_ptr.push_back(&*b);

  // compute self-intersections filtered out by boxes
  CGAL::internal::Intersect_facets<TM,GeomTraits,Box,OutputIterator> intersect_facets(tmesh, out, geom_traits);
  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}

/**
 * \ingroup PkgPolygonMeshProcessing
 * Checks if a polygon mesh is self-intersecting.
 * Depends on \ref PkgBoxIntersectionDSummary
 * @pre @a CGAL::is_pure_triangle(tmesh)
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam TriangleMesh a model of `FaceListGraph` (possibly a %CGAL polyhedron)
 *
 * @param tmesh TriangleMesh to be tested
 * @param geom_traits traits class providing intersection test primitives
 *
 * @return true if `tmesh` is self-intersecting
 */
template <class GeomTraits, class TriangleMesh>
bool is_self_intersecting(const TriangleMesh& tmesh,
                          const GeomTraits& geom_traits = GeomTraits())
{
  CGAL_precondition(CGAL::is_pure_triangle(tmesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    self_intersections<GeomTraits>(tmesh, OutputIterator(), geom_traits); 
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& ) 
  { return true; }

  return false;
}

}// end namespace Polygon_mesh_processing

}// namespace CGAL

#endif // CGAL_SELF_INTERSECTION_POLYHEDRON_3
