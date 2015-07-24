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
#include <boost/range.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

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
      Triangle t1 = triangle_functor( get(m_vpmap,target(h,m_tmesh)),
                                      get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                      get(m_vpmap, target(next(next(h,m_tmesh),m_tmesh),m_tmesh)));
      Triangle t2 = triangle_functor( get(m_vpmap, target(v,m_tmesh)),
                                      get(m_vpmap, target(next(v,m_tmesh),m_tmesh)),
                                      get(m_vpmap, target(next(next(v,m_tmesh),m_tmesh),m_tmesh)));
      
      Segment s1 = segment_functor( get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                    get(m_vpmap, target(next(next(h,m_tmesh),m_tmesh),m_tmesh)));
      Segment s2 = segment_functor( get(m_vpmap, target(next(v,m_tmesh),m_tmesh)),
                                    get(m_vpmap, target(next(next(v,m_tmesh),m_tmesh),m_tmesh)));
      
      if(do_intersect_3_functor(t1,s2)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      } else if(do_intersect_3_functor(t2,s1)){
        *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
      }
      return;
    }
    
    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap, target(h,m_tmesh)),
                                    get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                    get(m_vpmap, target(next(next(h,m_tmesh),m_tmesh),m_tmesh)));
    Triangle t2 = triangle_functor( get(m_vpmap, target(g,m_tmesh)),
                                    get(m_vpmap, target(next(g,m_tmesh),m_tmesh)),
                                    get(m_vpmap, target(next(next(g,m_tmesh),m_tmesh),m_tmesh)));
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
template <class TriangleMesh
        , class FaceRange
        , class OutputIterator
        , class NamedParameters
>
OutputIterator
self_intersections( const FaceRange& face_range,
                    const TriangleMesh& tmesh,
                    OutputIterator out,
                    const NamedParameters& np);
#endif

/** 
 * \ingroup PkgPolygonMeshProcessing
 * detects and records self-intersections of a triangulated surface mesh.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam TriangleMesh a model of `FaceListGraph` that has an internal property map
*         for `CGAL::vertex_point_t`
 * @tparam OutputIterator a model of `OutputIterator` holding objects of type 
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `SelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return `out`
 */
template <class TriangleMesh
        , class OutputIterator
#ifdef DOXYGEN_RUNNING
        , class NamedParameters
#else //avoid ambiguity with self_intersections(faces, tmesh, out)
        , class P, class T, class R
#endif
>
OutputIterator
self_intersections(const TriangleMesh& tmesh
                 , OutputIterator out
#ifdef DOXYGEN_RUNNING
                 , const NamedParameters& np)
#else
                 , const pmp_bgl_named_params<P,T,R>& np)
#endif
{
  return self_intersections(faces(tmesh), tmesh, out, np);
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh, class OutputIterator>
OutputIterator
self_intersections(const TriangleMesh& tmesh, OutputIterator out)
{
  return self_intersections(tmesh, out,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
/// \endcond

/*!
 * \ingroup PkgPolygonMeshProcessing
 * Same as above but the self-intersections reported
 * are only limited to the faces in `face_range`.
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
 *  model of `Range`.
 * Its iterator type is `RandomAccessIterator`.
 * @tparam TriangleMesh a model of `FaceListGraph` that has an internal property map
 *         for `CGAL::vertex_point_t`
 * @tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param face_range the range of faces to check for self-intersection.
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `SelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd

 */
template <class TriangleMesh
        , class FaceRange
        , class OutputIterator
        , class NamedParameters
>
OutputIterator
self_intersections( const FaceRange& face_range,
                    const TriangleMesh& tmesh,
                    OutputIterator out,
                    const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor> Box;

  // make one box per facet
  std::vector<Box> boxes;
  boxes.reserve(
    std::distance( boost::begin(face_range), boost::end(face_range) )
  );

  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_const_pmap(get_param(np, boost::vertex_point),
                                           tmesh,
                                           boost::vertex_point);

  BOOST_FOREACH(face_descriptor f, face_range)
  {
    boxes.push_back(Box( get(vpmap, target(halfedge(f,tmesh),tmesh)).bbox()
      + get(vpmap, target(next(halfedge(f, tmesh), tmesh), tmesh)).bbox()
      + get(vpmap, target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh)).bbox(),
    f));
  }
  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(num_faces(tmesh));

  BOOST_FOREACH(Box& b, boxes)
    box_ptr.push_back(&b);

  // compute self-intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;
  CGAL::internal::Intersect_facets<TM,GeomTraits,Box,OutputIterator,VertexPointMap>
    intersect_facets(tmesh, out, vpmap,
      boost::choose_param(get_param(np, geom_traits), GeomTraits()));

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh
        , class FaceRange
        , class OutputIterator
>
OutputIterator self_intersections(const FaceRange& face_range,
                                  const TriangleMesh& tmesh,
                                  OutputIterator out)
{
  return self_intersections(face_range, tmesh, out,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
/// \endcond

/**
 * \ingroup PkgPolygonMeshProcessing
 * tests if a triangulated surface mesh self-intersects.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam TriangleMesh a model of `FaceListGraph` that has an internal property map
 *         for `CGAL::vertex_point_t`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param tmesh the triangulated surface mesh to be tested
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh` \cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `SelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return true if `tmesh` self-intersects
 */
template <class TriangleMesh
        , class NamedParameters
          >
bool does_self_intersect(const TriangleMesh& tmesh
                        , const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    self_intersections(tmesh, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& ) 
  { return true; }

  return false;
}

/// \cond SKIP_IN_MANUAL
template <class TriangleMesh>
bool does_self_intersect(const TriangleMesh& tmesh)
{
  return does_self_intersect(tmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}
/// \endcond

}// end namespace Polygon_mesh_processing

}// namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
