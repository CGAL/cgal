// Copyright (c) 2008 INRIA Sophia-Antipolis (France).
// Copyright (c) 2008-2015 GeometryFactory (France).
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

#ifndef CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
#define CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#include <vector>
#include <exception>
#include <boost/foreach.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

namespace CGAL {
namespace internal {
template <class TM,//TriangleMesh
          class Kernel,
          class Box,
          class OutputIterator,
          class VertexPointMap>
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
  const VertexPointMap m_vpmap;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  
  Intersect_facets(const TM& tmesh, OutputIterator it, VertexPointMap vpmap, const Kernel& kernel)
    : 
    m_tmesh(tmesh),
    m_vpmap(vpmap),
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
      Triangle t1 = triangle_functor( m_vpmap[target(h,m_tmesh)],
                                      m_vpmap[target(next(h,m_tmesh),m_tmesh)],
                                      m_vpmap[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
      Triangle t2 = triangle_functor( m_vpmap[target(v,m_tmesh)],
                                      m_vpmap[target(next(v,m_tmesh),m_tmesh)],
                                      m_vpmap[target(next(next(v,m_tmesh),m_tmesh),m_tmesh)]);
      
      Segment s1 = segment_functor( m_vpmap[target(next(h,m_tmesh),m_tmesh)],
                                    m_vpmap[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
      Segment s2 = segment_functor( m_vpmap[target(next(v,m_tmesh),m_tmesh)],
                                    m_vpmap[target(next(next(v,m_tmesh),m_tmesh),m_tmesh)]);
      
      if(do_intersect_3_functor(t1,s2)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      } else if(do_intersect_3_functor(t2,s1)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      }
      return;
    }
    
    // check for geometric intersection
    Triangle t1 = triangle_functor( m_vpmap[target(h,m_tmesh)],
                                    m_vpmap[target(next(h,m_tmesh),m_tmesh)],
                                    m_vpmap[target(next(next(h,m_tmesh),m_tmesh),m_tmesh)]);
    Triangle t2 = triangle_functor( m_vpmap[target(g,m_tmesh)],
                                    m_vpmap[target(next(g,m_tmesh),m_tmesh)],
                                    m_vpmap[target(next(next(g,m_tmesh),m_tmesh),m_tmesh)]);
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

#ifndef DOXYGEN_RUNNING
template <class GeomTraits
        , class TriangleMesh
        , class FaceRange
        , class OutputIterator
        , class VertexPointMap
>
OutputIterator
self_intersections( const FaceRange& face_range,
                    const TriangleMesh& tmesh,
                    OutputIterator out,
                    VertexPointMap vpmap,
                    const GeomTraits& geom_traits);
#endif

/** 
 * \ingroup PkgPolygonMeshProcessing
 * detects and reports self-intersections of a triangle mesh.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_pure_triangle(tmesh)`
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam OutputIterator a model of `OutputIterator` holding objects of type 
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with
 `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
 `Kernel::Point_3` as value type.
 *
 * @param tmesh triangle mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect
 * @param vpmap the property map with the points associated to the vertices of `pmesh`
 * @param geom_traits geometric traits class providing intersection test primitives
 *
 * \todo code : named parameters
 *
 * @return `out`
 */
template <class GeomTraits
        , class TriangleMesh
        , class OutputIterator
        , class VertexPointMap
>
OutputIterator
self_intersections(const TriangleMesh& tmesh
                 , OutputIterator out
                 , VertexPointMap vpmap = VertexPointMap()
                 , const GeomTraits& geom_traits = GeomTraits())
{
  return self_intersections(faces(tmesh), tmesh, out, vpmap, geom_traits);
}

/*!
 * \ingroup PkgPolygonMeshProcessing
 * Same as the previous overload but the self-intersection reported
 * are only limited to the faces in `face_range`.
 *  \tparam FaceRange a range of face descriptors (nested type `iterator`,
 *  `const_iterator` and functions `begin()`, `end()` must be available.
 *  \param face_range the set of faces to check for self-intersection.
 * \todo code: is it a bad idea to use std::distance?
 */
template <class GeomTraits
        , class TriangleMesh
        , class FaceRange
        , class OutputIterator
        , class VertexPointMap
>
OutputIterator
self_intersections( const FaceRange& face_range,
                    const TriangleMesh& tmesh,
                    OutputIterator out,
                    VertexPointMap vpmap,
                    const GeomTraits& geom_traits)
{
  CGAL_precondition(CGAL::is_pure_triangle(tmesh));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor> Box;

  // make one box per facet
  std::vector<Box> boxes;
  boxes.reserve(
    std::distance(face_range.begin(), face_range.end())
  );

  BOOST_FOREACH(face_descriptor f, face_range)
  {
    boxes.push_back(Box( vpmap[target(halfedge(f,tmesh),tmesh)].bbox() +
                         vpmap[target(next(halfedge(f, tmesh), tmesh), tmesh)].bbox() +
                         vpmap[target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh)].bbox(),
    f));
  }
  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(num_faces(tmesh));

  BOOST_FOREACH(Box& b, boxes)
    box_ptr.push_back(&b);

  // compute self-intersections filtered out by boxes
  CGAL::internal::Intersect_facets<TM,GeomTraits,Box,OutputIterator,VertexPointMap>
    intersect_facets(tmesh, out, vpmap, geom_traits);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}

/**
 * \ingroup PkgPolygonMeshProcessing
 * checks if a triangle mesh is self-intersecting.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_pure_triangle(tmesh)`
 *
 * @tparam GeomTraits a model of `SelfIntersectionTraits`
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam VertexPointMap a model of `ReadablePropertyMap` with
 `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and
 `Kernel::Point_3` as value type.
 *
 * @param tmesh triangle mesh to be tested
 * @param geom_traits traits class providing intersection test primitives
 * @param vpmap the property map with the points associated to the vertices of `pmesh`
 *
 * @return true if `tmesh` is self-intersecting
 * \todo code : named parameters
 */
template <class GeomTraits
        , class TriangleMesh
        , class VertexPointMap
//         = typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type
          >
bool is_self_intersecting(const TriangleMesh& tmesh
                        , VertexPointMap vpmap = VertexPointMap()
                        , const GeomTraits& geom_traits = GeomTraits())
{
  CGAL_precondition(CGAL::is_pure_triangle(tmesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    self_intersections<GeomTraits>(tmesh, OutputIterator(), vpmap, geom_traits);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& ) 
  { return true; }

  return false;
}

}// end namespace Polygon_mesh_processing

}// namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
