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

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>
#include <CGAL/exceptions.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#endif

#include <boost/function_output_iterator.hpp>

#include <exception>
#include <sstream>
#include <type_traits>
#include <vector>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

// Checks for 'real' intersections, i.e. not simply a shared vertex or edge
template <class GT, class TM, class VPM>
bool do_faces_intersect(typename boost::graph_traits<TM>::halfedge_descriptor h,
                        typename boost::graph_traits<TM>::halfedge_descriptor g,
                        const TM& tmesh,
                        const VPM vpmap,
                        const typename GT::Construct_segment_3& construct_segment,
                        const typename GT::Construct_triangle_3& construct_triangle,
                        const typename GT::Do_intersect_3& do_intersect)
{
  typedef typename boost::graph_traits<TM>::vertex_descriptor                 vertex_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor               halfedge_descriptor;

  typedef typename GT::Segment_3                                              Segment;
  typedef typename GT::Triangle_3                                             Triangle;

  vertex_descriptor hv[3], gv[3];
  hv[0] = target(h, tmesh);
  hv[1] = target(next(h, tmesh), tmesh);
  hv[2] = source(h, tmesh);

  gv[0] = target(g, tmesh);
  gv[1] = target(next(g, tmesh), tmesh);
  gv[2] = source(g, tmesh);

  halfedge_descriptor opp_h;

  // check for shared edge
  for(unsigned int i=0; i<3; ++i)
  {
    opp_h = opposite(h, tmesh);
    if(face(opp_h, tmesh) == face(g, tmesh))
    {
      // there is an intersection if the four points are coplanar and the triangles overlap
      if(CGAL::coplanar(get(vpmap, hv[i]),
                        get(vpmap, hv[(i+1)%3]),
                        get(vpmap, hv[(i+2)%3]),
                        get(vpmap, target(next(opp_h, tmesh), tmesh))) &&
         CGAL::coplanar_orientation(get(vpmap, hv[(i+2)%3]),
                                    get(vpmap, hv[i]),
                                    get(vpmap, hv[(i+1)%3]),
                                    get(vpmap, target(next(opp_h, tmesh), tmesh)))
           == CGAL::POSITIVE)
      {
        return true;
      }
      else
      {
        // there is a shared edge but no intersection
        return false;
      }
    }

    h = next(h, tmesh);
  }

  // check for shared vertex --> maybe intersection, maybe not
  int i(0), j(0);
  bool shared = false;
  for(; i<3 && (! shared); ++i)
  {
    for(j=0; j<3 && (! shared); ++j)
    {
      if(hv[i] == gv[j])
      {
        shared = true;
        break;
      }
    }

    if(shared)
      break;
  }

  if(shared)
  {
    // found shared vertex:
    CGAL_assertion(hv[i] == gv[j]);

    // geometric check if the opposite segments intersect the triangles
    Triangle t1 = construct_triangle(get(vpmap, hv[0]), get(vpmap, hv[1]), get(vpmap, hv[2]));
    Triangle t2 = construct_triangle(get(vpmap, gv[0]), get(vpmap, gv[1]), get(vpmap, gv[2]));

    Segment s1 = construct_segment(get(vpmap, hv[(i+1)%3]), get(vpmap, hv[(i+2)%3]) );
    Segment s2 = construct_segment(get(vpmap, gv[(j+1)%3]), get(vpmap, gv[(j+2)%3]));

    if(do_intersect(t1, s2))
      return true;
    else if(do_intersect(t2, s1))
      return true;

    return false;
  }

  // check for geometric intersection
  Triangle t1 = construct_triangle(get(vpmap, hv[0]), get(vpmap, hv[1]), get(vpmap, hv[2]));
  Triangle t2 = construct_triangle(get(vpmap, gv[0]), get(vpmap, gv[1]), get(vpmap, gv[2]));
  if(do_intersect(t1, t2))
    return true;

  return false;
};

template <class Box, class TM, class VPM, class GT,
          class OutputIterator>
struct Strict_intersect_faces // meaning, not just a shared subface
{
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor;

  mutable OutputIterator m_iterator;
  const TM& m_tmesh;
  const VPM m_vpmap;
  typename GT::Construct_segment_3 m_construct_segment;
  typename GT::Construct_triangle_3 m_construct_triangle;
  typename GT::Do_intersect_3 m_do_intersect;

  Strict_intersect_faces(const TM& tmesh, VPM vpmap, const GT& gt, OutputIterator it)
    :
      m_iterator(it),
      m_tmesh(tmesh),
      m_vpmap(vpmap),
      m_construct_segment(gt.construct_segment_3_object()),
      m_construct_triangle(gt.construct_triangle_3_object()),
      m_do_intersect(gt.do_intersect_3_object())
  {}

  void operator()(const Box* b, const Box* c) const
  {
    halfedge_descriptor h = halfedge(b->info(), m_tmesh);
    halfedge_descriptor g = halfedge(c->info(), m_tmesh);

    if(do_faces_intersect<GT>(h, g, m_tmesh, m_vpmap, m_construct_segment, m_construct_triangle, m_do_intersect))
      *m_iterator++ = std::make_pair(b->info(), c->info());
  }
};

#ifdef CGAL_LINKED_WITH_TBB
// The functor doing all geometric tests in parallel
template <class FacePairs, class DoIntersectVector,
          class TM, class VPM, class GT>
struct Concurrent_face_intersection_tester
{
  typedef typename boost::graph_traits<TM>::halfedge_descriptor   halfedge_descriptor;
  typedef typename boost::graph_traits<TM>::face_descriptor       face_descriptor;

  const TM& m_tmesh;
  const VPM m_vpmap;
  const FacePairs& m_face_pairs;
  DoIntersectVector& m_do_intersect_vector;
  typename GT::Construct_segment_3 m_construct_segment;
  typename GT::Construct_triangle_3 m_construct_triangle;
  typename GT::Do_intersect_3 m_do_intersect;

  Concurrent_face_intersection_tester(const FacePairs& face_pairs,
                                      DoIntersectVector& do_intersect_vector,
                                      const TM& tmesh,
                                      VPM vpmap,
                                      const GT& gt)
    :
      m_tmesh(tmesh),
      m_vpmap(vpmap),
      m_face_pairs(face_pairs),
      m_do_intersect_vector(do_intersect_vector),
      m_construct_segment(gt.construct_segment_3_object()),
      m_construct_triangle(gt.construct_triangle_3_object()),
      m_do_intersect(gt.do_intersect_3_object())
  {}

  void operator()(const tbb::blocked_range<std::size_t> &r) const
  {
    for(std::size_t ri = r.begin(); ri != r.end(); ++ri)
      this->operator()(ri);
  }

  void operator()(std::size_t ri) const
  {
    const std::pair<face_descriptor,face_descriptor>& ff = m_face_pairs[ri];
    halfedge_descriptor h = halfedge(ff.first, m_tmesh), g = halfedge(ff.second, m_tmesh);

    if(do_faces_intersect<GT>(h, g, m_tmesh, m_vpmap, m_construct_segment, m_construct_triangle, m_do_intersect))
      m_do_intersect_vector[ri] = true;
  }
};

// This filter does not filter anything, but simply forwards intersecting pair information so that
// filtering can be done outside of the call to 'box_intersection_d'
template <class OutputIterator>
struct All_faces_filter
{
  All_faces_filter(OutputIterator it) :  m_iterator(it) { }

  template <class Box>
  void operator()(const Box* b, const Box* c) const { *m_iterator++ = std::make_pair(b->info(), c->info()); }

  mutable OutputIterator m_iterator;
};
#endif

template <class ConcurrencyTag,
          class TriangleMesh,
          class FaceRange,
          class FacePairOutputIterator,
          class NamedParameters>
FacePairOutputIterator
self_intersections_impl(const FaceRange& face_range,
                        const TriangleMesh& tmesh,
                        FacePairOutputIterator out,
                        const bool throw_on_SI,
                        const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  typedef TriangleMesh                                                                   TM;
  typedef typename boost::graph_traits<TM>::face_descriptor                              face_descriptor;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS                                  Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor, Box_policy> Box;

  typedef typename GetGeomTraits<TM, NamedParameters>::type                              GT;
  GT gt = choose_parameter(get_parameter(np, internal_np::geom_traits), GT());

  typedef typename GetVertexPointMap<TM, NamedParameters>::const_type                    VPM;
  VPM vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(boost::vertex_point, tmesh));

  const std::ptrdiff_t cutoff = 2000;

  // make one box per face
  std::vector<Box> boxes;
  boxes.reserve(std::distance(std::begin(face_range), std::end(face_range)));

  // This loop is very cheap, so there is hardly anything to gain from parallelizing it
  for(face_descriptor f : face_range)
  {
    typename boost::property_traits<VPM>::reference
      p = get(vpmap, target(halfedge(f,tmesh),tmesh)),
      q = get(vpmap, target(next(halfedge(f, tmesh), tmesh), tmesh)),
      r = get(vpmap, target(next(next(halfedge(f, tmesh), tmesh), tmesh), tmesh));

    // tiny fixme: if f is degenerate, we might still have a real intersection between f
    // and another face f', but right now we are not creating a box for f and thus not returning those
    if(collinear(p, q, r))
    {
      if(throw_on_SI)
        throw CGAL::internal::Throw_at_output_exception();
      else
        *out++= std::make_pair(f, f);
    }
    else
    {
      boxes.push_back(Box(p.bbox() + q.bbox() + r.bbox(), f));
    }
  }

  // generate box pointers
  std::vector<const Box*> box_ptr;
  box_ptr.reserve(boxes.size());

  for(Box& b : boxes)
    box_ptr.push_back(&b);

  // In case we are throwing, like in `does_self_intersect()`, we keep the geometric test to throw ASAP.
  // This is obviously not optimal if there are no or few self-intersections: it would be a greater speed-up
  // to do the same as for `self_intersections()`. However, doing like `self_intersections()` would
  // be a major slow-down over sequential code if there are a lot of self-intersections...
  typedef boost::function_output_iterator<CGAL::internal::Throw_at_output>               Throwing_output_iterator;
  typedef internal::Strict_intersect_faces<Box, TM, VPM, GT, Throwing_output_iterator>   Throwing_filter;
  Throwing_filter throwing_filter(tmesh, vpmap, gt, Throwing_output_iterator());

#if !defined(CGAL_LINKED_WITH_TBB)
  CGAL_static_assertion_msg (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    // (A) Parallel: Write all pairs of faces with intersecting bbox
    typedef tbb::concurrent_vector<std::pair<face_descriptor, face_descriptor> >  Face_pairs;
    typedef std::back_insert_iterator<Face_pairs>                                 Face_pairs_back_inserter;

    Face_pairs face_pairs;
    internal::All_faces_filter<Face_pairs_back_inserter> all_faces_filter(std::back_inserter(face_pairs));

    if(throw_on_SI)
      CGAL::box_self_intersection_d<ConcurrencyTag>(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
    else
      CGAL::box_self_intersection_d<ConcurrencyTag>(box_ptr.begin(), box_ptr.end(), all_faces_filter, cutoff);

    // (B) Parallel: Perform the geometric tests
    typedef std::vector<int>                                                      Do_intersect_vector;
    typedef internal::Concurrent_face_intersection_tester<Face_pairs, Do_intersect_vector, TM, VPM, GT> Tester;

    Do_intersect_vector do_intersect(face_pairs.size(), 0);
    Tester tester(face_pairs, do_intersect, tmesh, vpmap, gt);
    tbb::parallel_for(tbb::blocked_range<std::size_t>(0, face_pairs.size()), tester);

    // (C) Sequential: Copy from the concurent container to the output iterator
    for(std::size_t i=0; i<do_intersect.size(); ++i)
    {
      if(do_intersect[i])
        *out ++= face_pairs[i];
    }

    return out;
  }
#endif
  // Sequential version of the code
  // Compute self-intersections filtered out by boxes
  typedef internal::Strict_intersect_faces<Box, TM, VPM, GT, FacePairOutputIterator> Intersecting_faces_filter;
  Intersecting_faces_filter intersect_faces(tmesh, vpmap, gt, out);

  if(throw_on_SI)
    CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
  else
    CGAL::box_self_intersection_d(box_ptr.begin(), box_ptr.end(), intersect_faces, cutoff);

  return intersect_faces.m_iterator;
}

} // namespace internal

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
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam FaceRange a model of `ConstRange` with value type `boost::graph_traits<TriangleMesh>::%face_descriptor`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePairOutputIterator a model of `OutputIterator` holding objects of type
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
           class FacePairOutputIterator,
           class NamedParameters>
FacePairOutputIterator
self_intersections(const FaceRange& face_range,
                   const TriangleMesh& tmesh,
                         FacePairOutputIterator out,
                   const NamedParameters& np)
{
  return internal::self_intersections_impl<ConcurrencyTag>(face_range, tmesh, out, false /*don't throw*/, np);
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh,
          class FaceRange,
          class FacePairOutputIterator>
FacePairOutputIterator
self_intersections(const FaceRange& face_range,
                   const TriangleMesh& tmesh,
                         FacePairOutputIterator out)
{
  return self_intersections<ConcurrencyTag>(face_range, tmesh, out, CGAL::parameters::all_default());
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
 *                         Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePairOutputIterator a model of `OutputIterator` holding objects of type
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
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh,
          class FacePairOutputIterator,
          class CGAL_PMP_NP_TEMPLATE_PARAMETERS>
FacePairOutputIterator
self_intersections(const TriangleMesh& tmesh,
                         FacePairOutputIterator out,
                   const CGAL_PMP_NP_CLASS& np)
{
  return self_intersections<ConcurrencyTag>(faces(tmesh), tmesh, out, np);
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag, class TriangleMesh, class FacePairOutputIterator>
FacePairOutputIterator
self_intersections(const TriangleMesh& tmesh, FacePairOutputIterator out)
{
  return self_intersections<ConcurrencyTag>(faces(tmesh), tmesh, out, parameters::all_default());
}
/// \endcond

/**
 * \ingroup PMP_intersection_grp
 * tests if a set of faces of a triangulated surface mesh self-intersects.
 * This function depends on the package \ref PkgBoxIntersectionD
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
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
          class NamedParameters>
bool does_self_intersect(const FaceRange& face_range,
                         const TriangleMesh& tmesh,
                         const NamedParameters& np)
{
  CGAL_precondition(CGAL::is_triangle_mesh(tmesh));

  try
  {
    CGAL::Emptyset_iterator unused_out;
    internal::self_intersections_impl<ConcurrencyTag>(face_range, tmesh, unused_out, true /*throw*/, np);
  }
  catch(CGAL::internal::Throw_at_output_exception&)
  {
    return true;
  }

  return false;
}

/**
 * \ingroup PMP_intersection_grp
 * tests if a triangulated surface mesh self-intersects.
 * This function depends on the package \ref PkgBoxIntersectionD
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
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
          class TriangleMesh,
          class CGAL_PMP_NP_TEMPLATE_PARAMETERS>
bool does_self_intersect(const TriangleMesh& tmesh,
                         const CGAL_PMP_NP_CLASS& np)
{
  return does_self_intersect<ConcurrencyTag>(faces(tmesh), tmesh, np);
}

/// \cond SKIP_IN_MANUAL
template <class ConcurrencyTag = Sequential_tag, class TriangleMesh>
bool does_self_intersect(const TriangleMesh& tmesh)
{
  return does_self_intersect<ConcurrencyTag>(faces(tmesh), tmesh, CGAL::parameters::all_default());
}

template <class ConcurrencyTag = Sequential_tag, class FaceRange, class TriangleMesh>
bool does_self_intersect(const FaceRange& face_range,
                         const TriangleMesh& tmesh)
{
  return does_self_intersect<ConcurrencyTag>(face_range, tmesh, CGAL::parameters::all_default());
}
/// \endcond

}// namespace Polygon_mesh_processing
}// namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
