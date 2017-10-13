// Copyright (c) 2016 GeometryFactory (France).
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
// SPDX-License-Identifier: GPL-3.0+
//
//
// Author(s)     : Sebastien Loriot

#ifndef CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
#define CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H

#include <CGAL/license/Polygon_mesh_processing/corefinement.h>


#include <CGAL/Polygon_mesh_processing/internal/Corefinement/intersection_impl.h>
#include <boost/type_traits/is_same.hpp>

namespace CGAL {
namespace internal {

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_faces
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
  const TM& m_tmesh1;
  const VertexPointMap m_vpmap1;
  const TM& m_tmesh2;
  const VertexPointMap m_vpmap2;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_faces(const TM& mesh1, const TM& mesh2,
                   OutputIterator it,
                   VertexPointMap vpmap1, VertexPointMap vpmap2,
                   const Kernel& kernel)
    :
      m_tmesh1(mesh1),
      m_vpmap1(vpmap1),
      m_tmesh2(mesh2),
      m_vpmap2(vpmap2),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh1);
    halfedge_descriptor g = halfedge(c->info(), m_tmesh2);


    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap1, target(h,m_tmesh1)),
                                    get(m_vpmap1, target(next(h,m_tmesh1),m_tmesh1)),
                                    get(m_vpmap1, target(next(next(h,m_tmesh1),m_tmesh1),m_tmesh1)));

    Triangle t2 = triangle_functor( get(m_vpmap2, target(g,m_tmesh2)),
                                    get(m_vpmap2, target(next(g,m_tmesh2),m_tmesh2)),
                                    get(m_vpmap2, target(next(next(g,m_tmesh2),m_tmesh2),m_tmesh2)));
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_faces

template<class TM,
         class Kernel,
         class Box,
         class OutputIterator,
         class VertexPointMap>
struct Intersect_face_polyline
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
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;
  typedef typename boost::property_traits<Ppmap>::value_type Point;


  // members
  const TM& m_tmesh;
  const std::vector<face_descriptor>& faces;
  const VertexPointMap m_vpmap;
  const std::vector<Point>& polyline;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_face_polyline(const TM& mesh,
                           const std::vector<face_descriptor>& faces,
                           const std::vector<Point>& polyline,
                           OutputIterator it,
                           VertexPointMap vpmap,
                           const Kernel& kernel)
    :
      m_tmesh(mesh),
      faces(faces),
      m_vpmap(vpmap),
      polyline(polyline),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      segment_functor(kernel.construct_segment_3_object()),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(faces[b->info()], m_tmesh);


    // check for geometric intersection
    Triangle t = triangle_functor( get(m_vpmap, target(h,m_tmesh)),
                                   get(m_vpmap, target(next(h,m_tmesh),m_tmesh)),
                                   get(m_vpmap, target(next(next(h,m_tmesh),m_tmesh),m_tmesh)));

    Segment s = segment_functor(polyline[c->info()], polyline[c->info() + 1]);
    if(do_intersect_3_functor(t, s)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_face_polyline

template<class Polyline,
         class Kernel,
         class Box,
         class OutputIterator>
struct Intersect_polylines
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
  typedef typename Kernel::Point_3 Point;


  // members
  const std::vector<Point>& polyline1;
  const std::vector<Point>& polyline2;
  mutable OutputIterator  m_iterator;
  mutable bool            m_intersected;
  mutable boost::function_output_iterator<Output_iterator_with_bool> m_iterator_wrapper;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;




  Intersect_polylines(const Polyline& polyline1,
                      const Polyline& polyline2,
                      OutputIterator it,
                      const Kernel& kernel)
    :
      polyline1(polyline1),
      polyline2(polyline2),
      m_iterator(it),
      m_intersected(false),
      m_iterator_wrapper(Output_iterator_with_bool(&m_iterator, &m_intersected)),
      segment_functor(kernel.construct_segment_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const Box* b, const Box* c) const
  {


    // check for geometric intersection

    Segment s1 = segment_functor(polyline1[b->info()], polyline1[b->info() + 1]);
    Segment s2 = segment_functor(polyline2[c->info()], polyline2[c->info() + 1]);
    if(do_intersect_3_functor(s1, s2)){
      *m_iterator_wrapper++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_polylines

struct Throw_at_first_output {
  class Throw_at_first_output_exception: public std::exception
  { };

  template<class T>
  void operator()(const T& /* t */) const {
    throw Throw_at_first_output_exception();
  }
};
}// namespace internal

namespace Polygon_mesh_processing{

/*!
 * \ingroup PMP_intersection_grp
 * detects and records intersections between the specified
 * faces of two triangulated surface meshes.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(mesh1)`
 * \pre `CGAL::is_triangle_mesh(mesh2)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `Range`.
 * Its iterator type is `RandomAccessIterator`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param face_range1 the range of `mesh1` faces to check for intersections.
 * \param face_range2 the range of `mesh2` faces to check for intersections.
 * \param mesh1 the first triangulated surface mesh to be checked.
 * \param mesh2 the first triangulated surface mesh to be checked.
 * \param out output iterator to be filled with all pairs of faces that intersect
 * \param np1 optional sequence of \ref namedparameters for `mesh1`, among the ones listed below
 * \param np2 optional sequence of \ref namedparameters for `mesh2`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh
          , class FaceRange
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections( const FaceRange& face_range1,
               const FaceRange& face_range2,
               const TriangleMesh& mesh1,
               const TriangleMesh& mesh2,
               OutputIterator out,
               const NamedParameters& np1,
               const NamedParameters& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range1), boost::end(face_range1) )
        );
  boxes2.reserve(
        std::distance( boost::begin(face_range2), boost::end(face_range2) )
        );

  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap1 = boost::choose_param(get_param(np1, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, mesh1));
  VertexPointMap vpmap2 = boost::choose_param(get_param(np2, internal_np::vertex_point),
                                              get_const_property_map(boost::vertex_point, mesh2));

  BOOST_FOREACH(face_descriptor f, face_range1)
  {
    boxes1.push_back(Box( get(vpmap1, target(halfedge(f,mesh1),mesh1)).bbox()
                          + get(vpmap1, target(next(halfedge(f, mesh1), mesh1), mesh1)).bbox()
                          + get(vpmap1, target(next(next(halfedge(f, mesh1), mesh1), mesh1), mesh1)).bbox(),
                          f));
  }

  BOOST_FOREACH(face_descriptor f, face_range2)
  {
    boxes2.push_back(Box( get(vpmap2, target(halfedge(f,mesh2),mesh2)).bbox()
                          + get(vpmap2, target(next(halfedge(f, mesh2), mesh2), mesh2)).bbox()
                          + get(vpmap2, target(next(next(halfedge(f, mesh2), mesh2), mesh2), mesh2)).bbox(),
                          f));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(num_faces(mesh1));
  box2_ptr.reserve(num_faces(mesh2));

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

  CGAL::internal::Intersect_faces<TM,
      GeomTraits,
      Box,
      OutputIterator,
      VertexPointMap>
      Intersect_faces(mesh1, mesh2,
                       out,
                       vpmap1, vpmap2,
                       GeomTraits());

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           Intersect_faces,cutoff);
  return Intersect_faces.m_iterator;
}


/*!
 * \ingroup PMP_intersection_grp
 * detects and records intersections between the specified
 * faces of a triangulated surface mesh and a `Polyline`.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `Range`.
 * Its iterator type is `RandomAccessIterator`.
 * \tparam Polyline a `RandomAccessRange` of `CGAL::Point_3`. Those points must come from the same
 * geometric traits that templates`TriangleMesh`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param face_range the range of `mesh` faces to check for intersections.
 * \param polyline the `Polyline` to check for intersections.
 * \param mesh the triangulated surface mesh to be checked.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np optional sequence of \ref namedparameters for `mesh`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh
          , class FaceRange
          , class Polyline
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections( const FaceRange& face_range,
               const Polyline& polyline,
               const TriangleMesh& mesh,
               OutputIterator out,
               const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  typedef TriangleMesh TM;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;

  VertexPointMap vpmap = boost::choose_param(get_param(np, internal_np::vertex_point),
                                             get_const_property_map(boost::vertex_point, mesh));
  typedef typename boost::property_traits<VertexPointMap>::value_type Point;

  std::vector<face_descriptor> faces;
  faces.reserve(std::distance( boost::begin(face_range), boost::end(face_range) ));

  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;

  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(face_range), boost::end(face_range) )
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline), boost::end(polyline) ) - 1
        );


  BOOST_FOREACH(face_descriptor f, face_range)
  {
    faces.push_back(f);
    boxes1.push_back(Box( get(vpmap, target(halfedge(f,mesh),mesh)).bbox()
                          + get(vpmap, target(next(halfedge(f, mesh), mesh), mesh)).bbox()
                          + get(vpmap, target(next(next(halfedge(f, mesh), mesh), mesh), mesh)).bbox(),
                          faces.size()-1));
  }

  for(std::size_t i =0; i< polyline.size()-1; ++i)
  {
    Point p1 = polyline[i];
    Point p2 = polyline[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(boxes1.size());
  box2_ptr.reserve(boxes2.size());

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

  CGAL::internal::Intersect_face_polyline<TM,
      GeomTraits,
      Box,
      OutputIterator,
      VertexPointMap>
      Intersect_face_polyline(mesh,
                               faces,
                               polyline,
                               out,
                               vpmap,
                               GeomTraits());

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           Intersect_face_polyline, cutoff);
  return Intersect_face_polyline.m_iterator;
}

/*!
 * \ingroup PMP_intersection_grp
 * detects and records intersections between two `Polyline`s.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \tparam Polyline a `RandomAccessRange` of `CGAL::Point_3<Kernel>`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range.
 * \tparam Kernel a model of `Kernel`
 *
 * \param polyline1 the first polyline to check for intersections.
 * \param polyline2 the second polyline to check for intersections.
 * \param out output iterator to be filled with all pairs of segments that intersect
 * \param K an instance of `Kernel`
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template < class Polyline
           , class OutputIterator
           , class Kernel
           >
OutputIterator
intersections( const Polyline& polyline1,
               const Polyline& polyline2,
               OutputIterator out,
               const Kernel& K)
{
  typedef typename CGAL::Box_intersection_d::Box_with_info_d<double, 3, std::size_t> Box;
  typedef typename Kernel::Point_3 Point;
  // make one box per facet
  std::vector<Box> boxes1;
  std::vector<Box> boxes2;
  boxes1.reserve(
        std::distance( boost::begin(polyline1), boost::end(polyline1) ) - 1
        );

  boxes2.reserve(
        std::distance( boost::begin(polyline2), boost::end(polyline2) ) - 1
        );

  for(std::size_t i =0; i< polyline1.size()-1; ++i)
  {
    Point p1 = polyline1[i];
    Point p2 = polyline1[i+1];
    boxes1.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  for(std::size_t i =0; i< polyline2.size()-1; ++i)
  {
    Point p1 = polyline2[i];
    Point p2 = polyline2[i+1];
    boxes2.push_back(Box(p1.bbox() + p2.bbox(), i));
  }

  // generate box pointers
  std::vector<const Box*> box1_ptr;
  std::vector<const Box*> box2_ptr;
  box1_ptr.reserve(boxes1.size());
  box2_ptr.reserve(boxes2.size());

  BOOST_FOREACH(Box& b, boxes1)
      box1_ptr.push_back(&b);
  BOOST_FOREACH(Box& b, boxes2)
      box2_ptr.push_back(&b);

  // compute intersections filtered out by boxes

  CGAL::internal::Intersect_polylines<Polyline,
      Kernel,
      Box,
      OutputIterator>
      intersect_polylines(polyline1,
                          polyline2,
                          out,
                          K);

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_intersection_d(box1_ptr.begin(), box1_ptr.end(),
                           box2_ptr.begin(), box2_ptr.end(),
                           intersect_polylines, cutoff);
  return intersect_polylines.m_iterator;
}


/*!
 * \ingroup PMP_intersection_grp
 * detects and records intersections between two triangulated surface meshes.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * @pre `CGAL::is_triangle_mesh(mesh1)`
 * @pre `CGAL::is_triangle_mesh(mesh2)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param mesh1 the first triangulated surface mesh to be checked
 * \param mesh2 the second triangulated surface mesh to be checked
 * \param out output iterator to be filled with all pairs of faces that intersect
 * \param np1 optional sequence of \ref namedparameters for `mesh1`, among the ones listed below
 * \param np2 optional sequence of \ref namedparameters for `mesh2`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections(const TriangleMesh& mesh1,
              const TriangleMesh& mesh2,
              OutputIterator out,
              const NamedParameters& np1,
              const NamedParameters& np2)
{
  return intersections(faces(mesh1), faces(mesh2),
                       mesh1, mesh2, out, np1, np2);
}
/*!
 * \ingroup PMP_intersection_grp
 * detects and records intersections between a triangulated surface mesh
 * and a `Polyline`.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of `CGAL::Point_3`. Those points must come from the same
 * geometric traits that templates`TriangleMesh`.
 * \tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t, std::size_t>`. This OutputIterator will hold the position of the
 *  elements in their respective range.
 * \tparam NamedParameters a sequence of \ref namedparameters
 *
 * \param mesh the triangulated surface mesh to be checked.
 * \param polyline the `Polyline` to check for intersections.
 * \param out output iterator to be filled with all pairs of face-segment that intersect
 * \param np optional sequence of \ref namedparameters for `mesh`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 * \return `out`
 */
template <class TriangleMesh
          , class Polyline
          , class OutputIterator
          , class NamedParameters
          >
OutputIterator
intersections(const TriangleMesh& mesh,
              const Polyline& polyline,
              OutputIterator out,
              const NamedParameters& np)
{
  return intersections(faces(mesh), polyline, mesh, out, np);
}


/**
 * \ingroup PMP_intersection_grp
 * tests if two triangulated surface meshes intersect.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(mesh1)`
 * @pre `CGAL::is_triangle_mesh(mesh2)`
 *
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param mesh1 the first triangulated surface mesh to be tested
 * @param mesh2 the second triangulated surface mesh to be tested
 * @param np1 optional sequence of \ref namedparameters for `mesh1`, among the ones listed below
 * @param np2 optional sequence of \ref namedparameters for `mesh2`, among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return true if `mesh1` and `mesh2` intersect
 */
template <class TriangleMesh
          , class NamedParameters
          >
bool do_intersect(const TriangleMesh& mesh1
                  , const TriangleMesh& mesh2
                  , const NamedParameters& np1
                  , const NamedParameters& np2)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh1));
  CGAL_precondition(CGAL::is_triangle_mesh(mesh2));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    intersections(mesh1,mesh2, OutputIterator(), np1, np2);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 * tests if a triangulated surface mesh and
 * a `Polyline` intersect.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 * @pre `CGAL::is_triangle_mesh(mesh)`
 *
 * \tparam TriangleMesh a model of `FaceListGraph`
 * \tparam Polyline a `RandomAccessRange` of `CGAL::Point_3`. Those points must come from the same
 * geometric traits that templates`TriangleMesh`.
 * @tparam NamedParameters a sequence of \ref namedparameters
 *
 * @param mesh the triangulated surface mesh to be tested
 * @param polyline the `Polyline` to be tested.
 * @param np optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` should be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return true if `mesh` and `polyline` intersect
 */
template <class TriangleMesh
          , class Polyline
          , class NamedParameters
          >
bool do_intersect(const TriangleMesh& mesh
                  , const Polyline& polyline
                  , const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(mesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    intersections(mesh,polyline, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 * tests if two `Polyline`s intersect.
 * This function depends on the package \ref PkgBoxIntersectionDSummary
 *
 * \tparam Polyline a `RandomAccessRange` of `CGAL::Point_3<Kernel>`.
 * \tparam Kernel a model of `Kernel`.
 *
 * @param polyline1 the first `Polyline` to be tested.
 * @param polyline2 the second `Polyline` to be tested.
 * @param K an instance of `Kernel`.
 *
 * @return true if `polyline1` and `polyline2` intersect
 */
template < class Polyline
           , class Kernel
           >
bool do_intersect( const Polyline& polyline1
                  , const Polyline& polyline2,
                   const Kernel& K)
{
  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_first_output> OutputIterator;
    intersections(polyline1,polyline2, OutputIterator(), K);
  }
  catch( CGAL::internal::Throw_at_first_output::Throw_at_first_output_exception& )
  { return true; }

  return false;
}


/**
 * \ingroup PMP_corefinement_grp
 * computes the intersection of triangles of `tm1` and `tm2`. The output is a
 * set of polylines with all vertices but endpoints being of degree 2.
 *
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm1)` \endlink
 * \pre \link CGAL::Polygon_mesh_processing::does_self_intersect() `!CGAL::Polygon_mesh_processing::does_self_intersect(tm2)` \endlink
 *
 * @tparam TriangleMesh a model of `MutableFaceGraph`, `HalfedgeListGraph` and `FaceListGraph`
 * @tparam NamedParameters1 a sequence of \ref namedparameters
 * @tparam NamedParameters2 a sequence of \ref namedparameters
 * @tparam OutputIterator an output iterator in which `std::vector` of points
 *                        can be put. The point type is the one from the
 *                        vertex property map
 *
 * @param tm1 first input triangulated surface mesh
 * @param tm2 second input triangulated surface mesh
 * @param polyline_output output iterator of polylines. Each polyline will be
 *        given as a vector of points
 * @param throw_on_self_intersection if `true`, for each input triangle mesh,
 *        the set of triangles closed to the intersection of `tm1` and `tm2` will be
 *        checked for self-intersection and `CGAL::Corefinement::Self_intersection_exception`
 *        will be thrown if at least one is found.
 * @param np1 optional sequence of \ref namedparameters among the ones listed below
 * @param np2 optional sequence of \ref namedparameters among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamBegin{vertex_point_map}
 *    a property map with the points associated to the vertices of `tm1`
 *    (`tm2`). The two property map types must be the same.
 *    \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 */
template <class OutputIterator,
          class TriangleMesh,
          class NamedParameters1,
          class NamedParameters2 >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const NamedParameters1& np1,
                     const NamedParameters2& np2,
                     const bool throw_on_self_intersection=false)
{
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters1>::const_type Vpm;
  typedef typename GetVertexPointMap<TriangleMesh,
                                     NamedParameters2>::const_type Vpm2;
  CGAL_USE_TYPE(Vpm2);
  CGAL_assertion_code(
    static const bool same_vpm = (boost::is_same<Vpm,Vpm2>::value);)
  CGAL_static_assertion(same_vpm);

  Vpm vpm1 = choose_const_pmap(get_param(np1, internal_np::vertex_point),
                               tm1,
                               vertex_point);
  Vpm vpm2 = choose_const_pmap(get_param(np2, internal_np::vertex_point),
                               tm2,
                               vertex_point);

  Corefinement::Intersection_of_triangle_meshes<TriangleMesh,Vpm>
    functor(tm1, tm2, vpm1, vpm2);
  functor(polyline_output, throw_on_self_intersection, true);
  return polyline_output;
}

template <class OutputIterator,
          class TriangleMesh >
OutputIterator
surface_intersection(const TriangleMesh& tm1,
                     const TriangleMesh& tm2,
                     OutputIterator polyline_output,
                     const bool throw_on_self_intersection=false)
{
  return surface_intersection(tm1, tm2, polyline_output,
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    CGAL::Polygon_mesh_processing::parameters::all_default(),
    throw_on_self_intersection
  );
}

} } //end of namespace CGAL::Polygon_mesh_processing

#endif // CGAL_POLYGON_MESH_PROCESSING_INTERSECTION_H
