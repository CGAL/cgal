// Copyright (c) 2008 INRIA Sophia-Antipolis (France).
// Copyright (c) 2008-2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Pierre Alliez, Laurent Rineau, Ilker O. Yaz

// compute self-intersection of a CGAL triangle polyhedron mesh
// original code from Lutz Kettner

#ifndef CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
#define CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS

#include <CGAL/license/Polygon_mesh_processing/predicate.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/box_intersection_d.h>
#include <CGAL/intersections.h>
#include <CGAL/Bbox_3.h>

#ifdef CGAL_PMP_SI_DEBUG
#include <CGAL/Real_timer.h>
#endif

#include <CGAL/Kernel/global_functions_3.h>

#include <sstream>
#include <vector>
#include <exception>
#include <boost/range.hpp>

#include <boost/function_output_iterator.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/graph/graph_traits.hpp>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#endif

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {

namespace internal {

  template <class TM,//TriangleMesh
          class Kernel,
          class Box,
          class OutputIterator,
          class VertexPointMap>
struct Intersect_facets
{
// typedefs
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
    typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

// members
  const TM& m_tmesh;
  const VertexPointMap m_vpmap;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;

  Intersect_facets(const TM& tmesh, OutputIterator it, VertexPointMap vpmap, const Kernel& kernel)
    :
    m_tmesh(tmesh),
    m_vpmap(vpmap),
    m_iterator(it),
    segment_functor(kernel.construct_segment_3_object()),
    triangle_functor(kernel.construct_triangle_3_object()),
    do_intersect_3_functor(kernel.do_intersect_3_object())
  {}

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh);
    halfedge_descriptor g = halfedge(c->info(),m_tmesh);

    vertex_descriptor hv[3], gv[3];
    hv[0] = target(h, m_tmesh);
    hv[1] = target(next(h, m_tmesh), m_tmesh);
    hv[2] = source(h, m_tmesh);

    gv[0] = target(g, m_tmesh);
    gv[1] = target(next(g, m_tmesh), m_tmesh);
    gv[2] = source(g, m_tmesh);

    halfedge_descriptor opp_h;

    // check for shared egde
    for(unsigned int i=0; i<3; ++i){
      opp_h = opposite(h, m_tmesh);
      if(face(opp_h, m_tmesh) == c->info()){
        // there is an intersection if the four points are coplanar and
        // the triangles overlap
		  get(m_vpmap, hv[i]);
		  get(m_vpmap, hv[(i + 1) % 3]);
		  get(m_vpmap, hv[(i + 2) % 3]);
		  get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh));

        if(CGAL::coplanar(get(m_vpmap, hv[i]),
                          get(m_vpmap, hv[(i+1)%3]),
                          get(m_vpmap, hv[(i+2)%3]),
                          get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh))) &&
           CGAL::coplanar_orientation(get(m_vpmap, hv[(i+2)%3]),
                                      get(m_vpmap, hv[i]),
                                      get(m_vpmap, hv[(i+1)%3]),
                                      get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh)))
             == CGAL::POSITIVE){
          *m_iterator++ = std::make_pair(b->info(), c->info());
          return;
        } else { // there is a shared edge but no intersection
          return;
        }
      }
      h = next(h, m_tmesh);
    }

    // check for shared vertex --> maybe intersection, maybe not

    halfedge_descriptor v;

    int i(0),j(0);
    bool shared = false;
    for(; i < 3 && (! shared); ++i){
      for(j = 0; j < 3 && (! shared); ++j){
        if(hv[i]==gv[j]){
          shared = true;
          break;
        }
      }
      if (shared) {
        break;
      }
    }
    if(shared){
      // found shared vertex:
		assert(hv[i] == gv[j]);
      // geometric check if the opposite segments intersect the triangles
      Triangle t1 = triangle_functor( get(m_vpmap, hv[0]),
                                      get(m_vpmap, hv[1]),
                                      get(m_vpmap, hv[2]));
      Triangle t2 = triangle_functor( get(m_vpmap, gv[0]),
                                      get(m_vpmap, gv[1]),
                                      get(m_vpmap, gv[2]));

      Segment s1 = segment_functor( get(m_vpmap, hv[(i+1)%3]),
                                    get(m_vpmap, hv[(i+2)%3]) );
      Segment s2 = segment_functor( get(m_vpmap, gv[(j+1)%3]),
                                    get(m_vpmap, gv[(j+2)%3]));

      if(do_intersect_3_functor(t1,s2)){
        *m_iterator++ = std::make_pair(b->info(), c->info());
      } else if(do_intersect_3_functor(t2,s1)){
        *m_iterator++ = std::make_pair(b->info(), c->info());
      }
      return;
    }

    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap, hv[0]),
                                    get(m_vpmap, hv[1]),
                                    get(m_vpmap, hv[2]));
    Triangle t2 = triangle_functor( get(m_vpmap, gv[0]),
                                    get(m_vpmap, gv[1]),
                                    get(m_vpmap, gv[2]));
    if(do_intersect_3_functor(t1, t2)){
      *m_iterator++ = std::make_pair(b->info(), c->info());
    }
  } // end operator ()
}; // end struct Intersect_facets

#ifdef CGAL_LINKED_WITH_TBB
// The functor for testing only triangles that do not share an edge or vertex in parallel
// TO BE REMOVED?
template <class TM,//TriangleMesh
          class Kernel,
          class VertexPointMap>
struct TriangleTriangle {
  const TM& m_tmesh;
  const VertexPointMap m_vpmap;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;
  typedef typename Kernel::Triangle_3   Triangle;

  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  const std::vector<std::pair<face_descriptor,face_descriptor> >& seq_v_faces;
  std::vector<int>& dointersect;
  bool& result;

  TriangleTriangle(const std::vector<std::pair<face_descriptor,face_descriptor> >& seq_v_faces,
                   std::vector<int>& dointersect,
                   bool& result,
                   const TM& tmesh, VertexPointMap vpmap, const Kernel& kernel)
    : seq_v_faces(seq_v_faces),
      m_tmesh(tmesh),
      m_vpmap(vpmap),
      triangle_functor(kernel.construct_triangle_3_object()),
      do_intersect_3_functor(kernel.do_intersect_3_object()),
      dointersect(dointersect),
      result(result)
  {}

  void operator()(const tbb::blocked_range<std::size_t> &r) const
  {
    for (std::size_t i = r.begin(); i != r.end(); ++i) {
      const std::pair<face_descriptor,face_descriptor>& tt = seq_v_faces[i];

      halfedge_descriptor h = halfedge(tt.first, m_tmesh);
      halfedge_descriptor g = halfedge(tt.second ,m_tmesh);

      vertex_descriptor hv[3], gv[3];
      hv[0] = target(h, m_tmesh);
      hv[1] = target(next(h, m_tmesh), m_tmesh);
      hv[2] = source(h, m_tmesh);

      gv[0] = target(g, m_tmesh);
      gv[1] = target(next(g, m_tmesh), m_tmesh);
      gv[2] = source(g, m_tmesh);
      Triangle t1 = triangle_functor( get(m_vpmap, hv[0]),
                                      get(m_vpmap, hv[1]),
                                      get(m_vpmap, hv[2]));
      Triangle t2 = triangle_functor( get(m_vpmap, gv[0]),
                                      get(m_vpmap, gv[1]),
                                      get(m_vpmap, gv[2]));
      if(do_intersect_3_functor(t1, t2)){
        dointersect[i] = true;
        result = true;
      }
    }
  }
};

//TODO: use the same code for the linear pass?
// The functor for doing all geometric tests in parallel
template <class TM,//TriangleMesh
          class Kernel,
          class VertexPointMap>
struct AllPairs
{
// types
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;

// data members
  const TM& m_tmesh;
  const VertexPointMap m_vpmap;
  const std::vector<std::pair<face_descriptor,face_descriptor> >& seq_faces;
  std::vector<int>& dointersect;
  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;


  AllPairs(const std::vector<std::pair<face_descriptor,face_descriptor> >& seq_faces,
           std::vector<int>& dointersect,
           const TM& tmesh, VertexPointMap vpmap, const Kernel& kernel)
    : m_tmesh(tmesh)
    , m_vpmap(vpmap)
    , seq_faces(seq_faces)
    , dointersect(dointersect)
    , triangle_functor(kernel.construct_triangle_3_object())
    , do_intersect_3_functor(kernel.do_intersect_3_object())
  {}

  void operator()(const tbb::blocked_range<std::size_t> &r) const
  {
    for (std::size_t ri = r.begin(); ri != r.end(); ++ri) {

      this->operator()(ri);
    }
  }

  void operator()(std::size_t ri) const
  {
    const std::pair<face_descriptor,face_descriptor>& ff = seq_faces[ri];
    halfedge_descriptor h = halfedge(ff.first, m_tmesh), g = halfedge(ff.second, m_tmesh);

    vertex_descriptor hv[3], gv[3];
    hv[0] = target(h, m_tmesh);
    hv[1] = target(next(h, m_tmesh), m_tmesh);
    hv[2] = source(h, m_tmesh);

    gv[0] = target(g, m_tmesh);
    gv[1] = target(next(g, m_tmesh), m_tmesh);
    gv[2] = source(g, m_tmesh);

    halfedge_descriptor opp_h;

    // check for shared egde
    for(unsigned int i=0; i<3; ++i){
      opp_h = opposite(h, m_tmesh);
      if(face(opp_h, m_tmesh) == ff.second){
        // there is an intersection if the four points are coplanar and
        // the triangles overlap
		  get(m_vpmap, hv[i]);
		  get(m_vpmap, hv[(i + 1) % 3]);
		  get(m_vpmap, hv[(i + 2) % 3]);
		  get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh));

        if(CGAL::coplanar(get(m_vpmap, hv[i]),
                          get(m_vpmap, hv[(i+1)%3]),
                          get(m_vpmap, hv[(i+2)%3]),
                          get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh))) &&
           CGAL::coplanar_orientation(get(m_vpmap, hv[(i+2)%3]),
                                      get(m_vpmap, hv[i]),
                                      get(m_vpmap, hv[(i+1)%3]),
                                      get(m_vpmap, target(next(opp_h, m_tmesh), m_tmesh)))
             == CGAL::POSITIVE){
          dointersect[ri]=true;
          return;
        } else { // there is a shared edge but no intersection
          return;
        }
      }
      h = next(h, m_tmesh);
    }

    // check for shared vertex --> maybe intersection, maybe not

    halfedge_descriptor v;

    int i(0),j(0);
    bool shared = false;
    for(; i < 3 && (! shared); ++i){
      for(j = 0; j < 3 && (! shared); ++j){
        if(hv[i]==gv[j]){
          shared = true;
          break;
        }
      }
	  if (shared) {
		  break;
	  }
    }
    if(shared){
      // found shared vertex:
		assert(hv[i] == gv[j]);
      // geometric check if the opposite segments intersect the triangles
      Triangle t1 = triangle_functor( get(m_vpmap, hv[0]),
                                      get(m_vpmap, hv[1]),
                                      get(m_vpmap, hv[2]));
      Triangle t2 = triangle_functor( get(m_vpmap, gv[0]),
                                      get(m_vpmap, gv[1]),
                                      get(m_vpmap, gv[2]));

      Segment s1 = segment_functor( get(m_vpmap, hv[(i+1)%3]),
                                    get(m_vpmap, hv[(i+2)%3]) );
      Segment s2 = segment_functor( get(m_vpmap, gv[(j+1)%3]),
                                    get(m_vpmap, gv[(j+2)%3]));

      if(do_intersect_3_functor(t1,s2)){
        dointersect[ri] = true;
      } else if(do_intersect_3_functor(t2,s1)){
        dointersect[ri] = true;
      }
      return;
    }

    // check for geometric intersection
    Triangle t1 = triangle_functor( get(m_vpmap, hv[0]),
                                    get(m_vpmap, hv[1]),
                                    get(m_vpmap, hv[2]));
    Triangle t2 = triangle_functor( get(m_vpmap, gv[0]),
                                    get(m_vpmap, gv[1]),
                                    get(m_vpmap, gv[2]));
    if(do_intersect_3_functor(t1, t2)){
      dointersect[ri] = true;
    }
  }
};

// The functor for filtering pairs of faces that share an edge or vertex
template <class TM,//TriangleMesh
          class Box,
          class OutputIterator>
struct Incident_faces_filter
{
// typedefs
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;

// members
  const TM& m_tmesh;
  mutable OutputIterator  m_iterator;

  Incident_faces_filter(const TM& tmesh, OutputIterator it)
    : m_tmesh(tmesh),
      m_iterator(it)
  {}

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh);
    halfedge_descriptor g = halfedge(c->info(),m_tmesh);

    vertex_descriptor hv[3], gv[3];
    hv[0] = target(h, m_tmesh);
    hv[1] = target(next(h, m_tmesh), m_tmesh);
    hv[2] = source(h, m_tmesh);

    gv[0] = target(g, m_tmesh);
    gv[1] = target(next(g, m_tmesh), m_tmesh);
    gv[2] = source(g, m_tmesh);

    halfedge_descriptor opp_h;

    // ignore shared vertex as already dealt with

    halfedge_descriptor v;

    int i(0),j(0);
    for(; i < 3; ++i){
      for(j = 0; j < 3; ++j){
        if(hv[i]==gv[j]){
          return;
        }
      }
    }
    *m_iterator++ = std::make_pair(b->info(), c->info());

  } // end operator ()
};

// The functor that filters nothing
template <class OutputIterator>
struct All_faces_filter
{
  mutable OutputIterator  m_iterator;


  All_faces_filter(OutputIterator it)
    :  m_iterator(it)
  {}

  template <class Box>
  void operator()(const Box* b, const Box* c) const
  {
    *m_iterator ++ = std::make_pair(b->info(), c->info());

  } // end operator ()
};


// The functor that computes intersections for faces incident to a vertex
template <class TM,//TriangleMesh
          class Kernel,
          class Box,
          class OutputIterator,
          class VertexPointMap>
struct Intersect_facets_incident_to_vertex
{
// typedefs
  typedef typename Kernel::Point_3      Point;
  typedef typename Kernel::Segment_3    Segment;
  typedef typename Kernel::Triangle_3   Triangle;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::property_map<TM, boost::vertex_point_t>::const_type Ppmap;

// members
  const TM& m_tmesh;
  const VertexPointMap m_vpmap;
  mutable OutputIterator  m_iterator;

  typename Kernel::Construct_segment_3  segment_functor;
  typename Kernel::Construct_triangle_3 triangle_functor;
  typename Kernel::Do_intersect_3       do_intersect_3_functor;


  Intersect_facets_incident_to_vertex(const TM& tmesh, OutputIterator it, VertexPointMap vpmap, const Kernel& kernel)
    :
    m_tmesh(tmesh),
    m_vpmap(vpmap),
    m_iterator(it),
    segment_functor(kernel.construct_segment_3_object()),
    triangle_functor(kernel.construct_triangle_3_object()),
    do_intersect_3_functor(kernel.do_intersect_3_object())
  { }

  void operator()(const tbb::blocked_range<std::size_t> &r) const
  {
    for (std::size_t i = r.begin(); i != r.end(); ++i) {
      vertex_descriptor vd(i);

      halfedge_descriptor hd = halfedge(vd, m_tmesh), done(hd);
      if(is_border(hd, m_tmesh)){
          hd = prev(opposite(hd,m_tmesh),m_tmesh);
        }

      do { //for each hd with vd as target look at all other faces
        Point p = get(m_vpmap, vd);
        Point q = get(m_vpmap, target(next(hd,m_tmesh), m_tmesh));
        Point r = get(m_vpmap, source(hd,m_tmesh));

        Triangle th = triangle_functor(p, q, r);

        // first look at the directly shared edges
        halfedge_descriptor ohd = opposite(hd, m_tmesh);
        if((! is_border(ohd,m_tmesh)) && (hd < ohd)){
          vertex_descriptor ov = target(next(ohd,m_tmesh), m_tmesh);
          Point s = get(m_vpmap, ov);
          if(CGAL::coplanar(p,q,r,s) && CGAL::coplanar_orientation(r,p,q,s) == POSITIVE){
            *m_iterator++ = std::make_pair(face(hd,m_tmesh), face(ohd, m_tmesh));
          }
        }
        ohd = opposite(next(hd,m_tmesh), m_tmesh);
        if((! is_border(ohd,m_tmesh)) && (hd < ohd)){
          vertex_descriptor ov = target(next(ohd,m_tmesh), m_tmesh);
          Point s = get(m_vpmap, ov);
          if(CGAL::coplanar(p,q,r,s) && CGAL::coplanar_orientation(p,q,r,s) == POSITIVE){
            *m_iterator++ = std::make_pair(face(hd,m_tmesh), face(ohd, m_tmesh));
          }
        }
        halfedge_descriptor po = prev(opposite(hd,m_tmesh),m_tmesh);
        halfedge_descriptor op = opposite(next(hd,m_tmesh),m_tmesh);

        for(halfedge_descriptor start : halfedges_around_target(hd,m_tmesh)){
          if((start == hd) || (start == op) || (start == po) || is_border(start, m_tmesh) || (face(hd, m_tmesh) < face(start, m_tmesh))){
            continue;
          }
          Segment ss = segment_functor(get(m_vpmap, source(start,m_tmesh)),
                                       get(m_vpmap, target(next(start,m_tmesh), m_tmesh)));
          if(do_intersect_3_functor(th, ss)){
            *m_iterator++ = std::make_pair(face(hd,m_tmesh), face(start, m_tmesh));
          }
        }
        hd = prev(opposite(hd,m_tmesh),m_tmesh);
      }while(hd != done);
    }
  }
}; // end struct Intersect_facets_incident_to_vertex
#endif

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

/*!
 * \ingroup PMP_intersection_grp
 * collects intersections between a subset of faces of a triangulated surface mesh.
 * Two faces are said to intersect if the corresponding triangles intersect
 * and the intersection is not an edge nor a vertex incident to both faces.
 *
 * This function depends on the package \ref PkgBoxIntersectionD
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag` and `Parallel_tag`.
 * @tparam FaceRange range of `boost::graph_traits<TriangleMesh>::%face_descriptor`,
 *  model of `Range`.
 * Its iterator type is `RandomAccessIterator`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param face_range the range of faces to check for self-intersection.
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` must be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 */
template < class ConcurrencyTag = Sequential_tag,
           class TriangleMesh,
           class FaceRange,
           class OutputIterator,
           class NamedParameters>
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

  // make one box per face
  std::vector<Box> boxes;
  boxes.reserve(
    std::distance( boost::begin(face_range), boost::end(face_range) )
  );

  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                      get_const_property_map(boost::vertex_point, tmesh));

  for(face_descriptor f : face_range)
  {
    typename boost::property_traits<VertexPointMap>::reference
      p = get(vpmap, target(halfedge(f,tmesh),tmesh)),
      q = get(vpmap, target(next(halfedge(f, tmesh), tmesh), tmesh)),
      r = get(vpmap, target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh));

    if ( collinear(p, q, r) )
      *out++= std::make_pair(f,f);
    else
      boxes.push_back(Box(p.bbox() + q.bbox() + r.bbox(), f));
  }
  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(boxes.size());

  for(Box& b : boxes)
    box_ptr.push_back(&b);

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(boost::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if (boost::is_convertible<ConcurrencyTag,Parallel_tag>::value)
  {
  // (Parallel version of the code)
    // (A) Sequentially write all pairs of faces with intersecting bbox into a std::vector
    std::ptrdiff_t cutoff = 2000;
    typedef std::vector<std::pair<face_descriptor,face_descriptor> >SeqV;
    typedef std::back_insert_iterator<SeqV> SeqVI;
    SeqV face_pairs;
    internal::All_faces_filter<SeqVI> all_faces_filter(std::back_inserter(face_pairs));
    CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),all_faces_filter,cutoff);

    // (B) Parallel: perform the geometric tests
    typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

    std::vector<int> dointersect(face_pairs.size(),0);

    CGAL::internal::AllPairs<TM,GeomTraits,VertexPointMap>
      all_pairs(face_pairs,
                dointersect,
                tmesh,
                vpmap,
                parameters::choose_parameter(parameters::get_parameter(np, internal_np::geom_traits), GeomTraits()));

    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, face_pairs.size()), all_pairs);

    // (C) Sequentially: Copy from the concurent container to the output iterator
    for(int i=0; i < dointersect.size(); ++i){
      if(dointersect[i])
        *out ++= face_pairs[i];
    }
    return out;
  }
#endif

  // (Sequential version of the code) compute self-intersections filtered out by boxes
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;
  CGAL::internal::Intersect_facets<TM,GeomTraits,Box,OutputIterator,VertexPointMap>
    intersect_facets(tmesh, out, vpmap,
      parameters::choose_parameter(parameters::get_parameter(np, internal_np::geom_traits), GeomTraits()));

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),intersect_facets,cutoff);
  return intersect_facets.m_iterator;
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh,
          class FaceRange,
          class OutputIterator>
OutputIterator
self_intersections(const FaceRange& face_range,
                   const TriangleMesh& tmesh,
                         OutputIterator out)
{
  return self_intersections(face_range, tmesh, out,
                            CGAL::Polygon_mesh_processing::parameters::all_default(),
                            ConcurrencyTag());
}
/// \endcond


/**
 * \ingroup PMP_intersection_grp
 * collects intersections between all the faces of a triangulated surface mesh.
 * Two faces are said to intersect if the corresponding triangles intersect
 * and the intersection is not an edge nor a vertex incident to both faces.
 *
 * This function depends on the package \ref PkgBoxIntersectionD
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag` and `Parallel_tag`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam OutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect.
              In case `tmesh` contains some degenerate faces, for each degenerate face `f` a pair `(f,f)`
              will be put in `out` before any other self intersection between non-degenerate faces.
              These are the only pairs where degenerate faces will be reported.
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` must be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return `out`
 */
template <class ConcurrencyTag = Sequential_tag
        , class TriangleMesh
        , class OutputIterator
        , class CGAL_PMP_NP_TEMPLATE_PARAMETERS>
OutputIterator
self_intersections(const TriangleMesh& tmesh,
                          OutputIterator out,
                   const CGAL_PMP_NP_CLASS& np)
{
  return self_intersections(faces(tmesh), tmesh, out, np);
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag, class TriangleMesh, class OutputIterator>
OutputIterator
self_intersections(const TriangleMesh& tmesh, OutputIterator out)
{
  return self_intersections<ConcurrencyTag>(faces(tmesh), tmesh, out, parameters::all_default());
}
/// \endcond

// experimental code for parallel processing on vertices
// (slower than general parallel on boxes). To be removed.
#if 0
template <class TriangleMesh
        , class FaceRange
        , class OutputIterator
        , class NamedParameters
>
OutputIterator
self_intersections( const FaceRange& face_range,
                    const TriangleMesh& tmesh,
                    OutputIterator out,
                    const NamedParameters& np,
                    int)
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
  VertexPointMap vpmap = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                                      get_const_property_map(boost::vertex_point, tmesh));

  for(face_descriptor f : face_range)
  {
    typename boost::property_traits<VertexPointMap>::reference
      p = get(vpmap, target(halfedge(f,tmesh),tmesh)),
      q = get(vpmap, target(next(halfedge(f, tmesh), tmesh), tmesh)),
      r = get(vpmap, target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh));

    if ( collinear(p, q, r) )
      *out++= std::make_pair(f,f);
    else
      boxes.push_back(Box(p.bbox() + q.bbox() + r.bbox(), f));
  }
  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(boxes.size());

  for(Box& b : boxes)
    box_ptr.push_back(&b);


  // (A) Parallel: look at each vertex and its incident faces
  //     Intersections are written into a concurrent vector
  typedef typename GetGeomTraits<TM, NamedParameters>::type GeomTraits;

  typedef tbb::concurrent_vector<std::pair<face_descriptor,face_descriptor> > CV;
  CV cv_faces;
  typedef std::back_insert_iterator<CV> CVI;

  CGAL::internal::Intersect_facets_incident_to_vertex<TM,GeomTraits,Box,CVI,VertexPointMap>
    intersect_facets_incident_to_vertex(tmesh,
                                        std::back_inserter(cv_faces), // result
                                        vpmap,
                                        parameters::choose_parameter(parameters::get_parameter(np, internal_np::geom_traits), GeomTraits()));
#ifdef CGAL_PMP_SI_DEBUG
  Real_timer rt;
  rt.start();
#endif

  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, num_vertices(tmesh)), intersect_facets_incident_to_vertex);

  // Copy from the concurent container to the output iterator
  for(typename CV::iterator it = cv_faces.begin(); it != cv_faces.end(); ++it){
    *out ++= *it;
  }
#ifdef CGAL_PMP_SI_DEBUG
  std::cout << "(A) Parallel: faces incident to each vertex " << rt.time() << "sec."<< std::endl;
  rt.reset();
#endif

  // (B) Sequential: Look at pairs of triangles with intersecting bbox
  //     Copy the pairs which do not share an edge or a vertex into a std::vector
  typedef std::vector<std::pair<face_descriptor,face_descriptor> >SeqV;
  SeqV seq_v_faces;
  typedef std::back_insert_iterator<SeqV> SeqVI;

  CGAL::internal::Incident_faces_filter<TM,Box,SeqVI>
                      incident_faces_filter(tmesh,
                                            std::back_inserter(seq_v_faces));

  std::ptrdiff_t cutoff = 2000;
  CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(),incident_faces_filter,cutoff);
#ifdef CGAL_PMP_SI_DEBUG
  std::cout << "(B) Sequential: boxes that intersect " << rt.time() << "sec."<< std::endl;
  rt.reset();
#endif

  // (C) Parallel: Perform geometric test on pairs not sharing an edge or a vertex
  std::vector<int> dointersect(seq_v_faces.size());
  bool result(false);
  CGAL::internal::TriangleTriangle<TM,GeomTraits,VertexPointMap> tt(seq_v_faces,
                                                                    dointersect,
                                                                    result,
                                                                    tmesh,
                                                                    vpmap,
                                                                    parameters::choose_parameter(parameters::get_parameter(np, internal_np::geom_traits), GeomTraits()));
  tbb::parallel_for(tbb::blocked_range<std::size_t>(0, seq_v_faces.size()), tt);
  if(result){
    for(int i = 0; i < seq_v_faces.size(); ++i){
      if(dointersect[i]){
        *out ++= seq_v_faces[i];
      }
    }
  }
#ifdef CGAL_PMP_SI_DEBUG
  std::cout << "(C) Parallel: Filter triangles that intersect " << rt.time() << "sec."<< std::endl;
#endif
  return out;
}
#endif // CGAL_LINKED_WITH_TBB

/**
 * \ingroup PMP_intersection_grp
 * tests if a triangulated surface mesh self-intersects.
 * This function depends on the package \ref PkgBoxIntersectionD
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag` and `Parallel_tag`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param tmesh the triangulated surface mesh to be tested
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` must be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `PMPSelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return `true` if `tmesh` self-intersects
 */
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh
        , class CGAL_PMP_NP_TEMPLATE_PARAMETERS
          >
bool does_self_intersect(const TriangleMesh& tmesh
                        , const CGAL_PMP_NP_CLASS& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    self_intersections<ConcurrencyTag>(tmesh, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& )
  { return true; }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 * tests if a set of faces of a triangulated surface mesh self-intersects.
 * This function depends on the package \ref PkgBoxIntersectionD
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam FaceRange a range of `face_descriptor`
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
 *
 * @param face_range the set of faces to test for self-intersection
 * @param tmesh the triangulated surface mesh to be tested
 * @param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *    \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `tmesh`.
 *   If this parameter is omitted, an internal property map for
 *   `CGAL::vertex_point_t` must be available in `TriangleMesh`\cgalParamEnd
 *    \cgalParamBegin{geom_traits} an instance of a geometric traits class, model of `SelfIntersectionTraits` \cgalParamEnd
 * \cgalNamedParamsEnd
 *
 * @return `true` if the faces in `face_range` self-intersect
 */
template <class ConcurrencyTag = Sequential_tag,
          class FaceRange,
          class TriangleMesh,
          class NamedParameters
          >
bool does_self_intersect(const FaceRange& face_range,
                         const TriangleMesh& tmesh,
                         const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  try
  {
    typedef boost::function_output_iterator<CGAL::internal::Throw_at_output> OutputIterator;
    self_intersections<ConcurrencyTag>(face_range, tmesh, OutputIterator(), np);
  }
  catch( CGAL::internal::Throw_at_output::Throw_at_output_exception& )
  { return true; }

  return false;
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag, class TriangleMesh>
bool does_self_intersect(const TriangleMesh& tmesh)
{
  return does_self_intersect<ConcurrencyTag>(tmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

template <class ConcurrencyTag = Sequential_tag, class FaceRange, class TriangleMesh>
bool does_self_intersect(const FaceRange& face_range,
                         const TriangleMesh& tmesh)
{
  return does_self_intersect<ConcurrencyTag>(face_range, tmesh,
    CGAL::Polygon_mesh_processing::parameters::all_default());
}

/// \endcond

}// end namespace Polygon_mesh_processing

}// namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
