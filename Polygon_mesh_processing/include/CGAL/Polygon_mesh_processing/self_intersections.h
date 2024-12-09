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

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/algorithm.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/exceptions.h>
#include <CGAL/intersections.h>
#include <CGAL/iterator.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/Random.h>
#include <CGAL/use.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#endif

#include <boost/iterator/function_output_iterator.hpp>
#include <boost/range/irange.hpp>

#include <exception>
#include <sstream>
#include <type_traits>
#include <typeinfo>
#include <vector>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <class TM>
struct Triangle_mesh_and_triangle_soup_wrapper
{
  typedef typename boost::graph_traits<TM>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TM>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TM>::halfedge_descriptor halfedge_descriptor; // private

  static void get_face_vertices(face_descriptor fd, std::array<vertex_descriptor,3>& vh, const TM& tm)
  {
    CGAL_assertion(boost::graph_traits<TM>::null_face() != fd);
    halfedge_descriptor h = halfedge(fd, tm);
    vh[0]=source(h, tm);
    vh[1]=target(h, tm);
    vh[2]=target(next(h, tm), tm);
  }

  static bool faces_have_a_shared_edge(face_descriptor f, face_descriptor g, std::array<vertex_descriptor, 4>& vh, const TM& tm)
  {
    CGAL_assertion(boost::graph_traits<TM>::null_face() != f);
    CGAL_assertion(boost::graph_traits<TM>::null_face() != g);
    halfedge_descriptor h=halfedge(f, tm);
    for(unsigned int i=0; i<3; ++i)
    {
      halfedge_descriptor opp_h = opposite(h, tm);
      if(face(opp_h, tm) == g)
      {
        vh[0]=source(h, tm);
        vh[1]=target(h, tm);
        vh[2]=target(next(h, tm), tm);
        vh[3]=target(next(opp_h, tm), tm);
        return true;
      }
      h = next(h, tm);
    }
    return false;
  }

  static bool is_pure_triangle(const TM& tm)
  {
    return is_triangle_mesh(tm);
  }
};

template <class PointRange, class TriangleRange>
struct Triangle_mesh_and_triangle_soup_wrapper<
    std::pair<const PointRange&,
              const TriangleRange&>>
{
  typedef std::size_t face_descriptor;
  typedef std::size_t vertex_descriptor;

  typedef std::pair<const PointRange&, const TriangleRange& > Soup;

  static void get_face_vertices(face_descriptor fd, std::array<vertex_descriptor,3>& vh, const Soup& soup)
  {
    const auto& face = soup.second[fd];
    vh[0]=face[0];
    vh[1]=face[1];
    vh[2]=face[2];
  }

  static bool faces_have_a_shared_edge(face_descriptor fd, face_descriptor gd, std::array<vertex_descriptor, 4>& vh, const Soup& soup)
  {
    const auto& f = soup.second[fd];
    const auto& g = soup.second[gd];

    for(unsigned int i=0; i<2; ++i) // no need to check f[2] if neither f[0] nor f[1] are shared
    {
      for(unsigned int j=0; j<3; ++j)
      {
        if (f[i]==g[j])
        {
          vh[0]=f[i];
          vh[1]=f[i+1];
          vh[2]=f[(i+2)%3];

          if (vh[1]==g[(j+1)%3])
          {
            vh[3]=g[(j+2)%3];
            return true;
          }
          if (vh[1]==g[(j+2)%3])
          {
            vh[3]=g[(j+1)%3];
            return true;
          }

          if (i==0)
          {
            vh[1]=f[i];
            vh[2]=f[(i+1)%3];
            vh[0]=f[(i+2)%3];
            if (vh[0]==g[(j+1)%3])
            {
              vh[3]=g[(j+2)%3];
              return true;
            }
            if (vh[0]==g[(j+2)%3])
            {
              vh[3]=g[(j+1)%3];
              return true;
            }
          }

          return false;
        }
      }
    }

    return false;
  }

  static bool is_pure_triangle(const Soup& soup)
  {
    for (const typename std::iterator_traits<typename TriangleRange::const_iterator>::value_type& t : soup.second)
      if (t.size()!=3)
        return false;
    return true;
  }
};

template<typename Output_iterator>
struct Throw_at_count_reached_functor {

  std::atomic<unsigned int>& counter;
  const unsigned int& maxval;

  Output_iterator out;

  Throw_at_count_reached_functor(std::atomic<unsigned int>& counter,
                                 const unsigned int& maxval,
                                 Output_iterator out)
    : counter(counter), maxval(maxval), out(out)
  {}

  template<class T>
  void operator()(const T& t )
  {
    *out++ = t;
    ++counter;
    if(counter >= maxval)
    {
      throw CGAL::internal::Throw_at_output_exception();
    }
  }
};

// Checks for 'real' intersections, i.e. not simply a shared vertex or edge
template <class GT, class TM, class VPM>
bool do_faces_intersect(typename Triangle_mesh_and_triangle_soup_wrapper<TM>::face_descriptor fh,
                        typename Triangle_mesh_and_triangle_soup_wrapper<TM>::face_descriptor fg,
                        const TM& tmesh,
                        const VPM vpmap,
                        const typename GT::Construct_segment_3& construct_segment,
                        const typename GT::Construct_triangle_3& construct_triangle,
                        const typename GT::Do_intersect_3& do_intersect)
{
  typedef Triangle_mesh_and_triangle_soup_wrapper<TM>       Wrapper;
  typedef typename Wrapper::vertex_descriptor     vertex_descriptor;

  typedef typename GT::Segment_3                            Segment;
  typedef typename GT::Triangle_3                          Triangle;


  std::array<vertex_descriptor, 3> hv, gv;
  Wrapper::get_face_vertices(fh, hv, tmesh);
  Wrapper::get_face_vertices(fg, gv, tmesh);

  // check for shared edge
  std::array<vertex_descriptor, 4> verts;
  if (Wrapper::faces_have_a_shared_edge(fh, fg, verts, tmesh))
  {
    if (verts[2]==verts[3]) return false; // only for a soup of triangles

    // there is an intersection if the four points are coplanar and the triangles overlap
    if(CGAL::coplanar(get(vpmap, verts[0]),
                      get(vpmap, verts[1]),
                      get(vpmap, verts[2]),
                      get(vpmap, verts[3])) &&
       CGAL::coplanar_orientation(get(vpmap, verts[0]),
                                  get(vpmap, verts[1]),
                                  get(vpmap, verts[2]),
                                  get(vpmap, verts[3]))
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
    const Triangle t1 = construct_triangle(get(vpmap, hv[0]), get(vpmap, hv[1]), get(vpmap, hv[2]));
    const Triangle t2 = construct_triangle(get(vpmap, gv[0]), get(vpmap, gv[1]), get(vpmap, gv[2]));

    const Segment s1 = construct_segment(get(vpmap, hv[(i+1)%3]), get(vpmap, hv[(i+2)%3]));
    const Segment s2 = construct_segment(get(vpmap, gv[(j+1)%3]), get(vpmap, gv[(j+2)%3]));

    if(do_intersect(t1, s2))
      return true;
    else if(do_intersect(t2, s1))
      return true;

    return false;
  }

  // check for geometric intersection
  const Triangle th = construct_triangle(get(vpmap, hv[0]), get(vpmap, hv[1]), get(vpmap, hv[2]));
  const Triangle tg = construct_triangle(get(vpmap, gv[0]), get(vpmap, gv[1]), get(vpmap, gv[2]));
  if(do_intersect(th, tg))
    return true;

  return false;
}

template <class Box, class TM, class VPM, class GT,
          class OutputIterator>
struct Strict_intersect_faces // "strict" as in "not sharing a subface"
{
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
    if(do_faces_intersect<GT>(b->info(), c->info(), m_tmesh, m_vpmap, m_construct_segment, m_construct_triangle, m_do_intersect))
      *m_iterator++ = std::make_pair(b->info(), c->info());
  }
};

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
  typedef TriangleMesh                                                                   TM;
  typedef Triangle_mesh_and_triangle_soup_wrapper<TM>                                    Wrapper;

  CGAL_precondition(Wrapper::is_pure_triangle(tmesh));

  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;
  using CGAL::parameters::is_default_parameter;

  typedef typename Wrapper::face_descriptor                                              face_descriptor;
  typedef typename Wrapper::vertex_descriptor                                            vertex_descriptor;

  typedef CGAL::Box_intersection_d::ID_FROM_BOX_ADDRESS                                  Box_policy;
  typedef CGAL::Box_intersection_d::Box_with_info_d<double, 3, face_descriptor, Box_policy> Box;

  typedef typename GetGeomTraits<TM, NamedParameters>::type                              GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  typedef GetVertexPointMap<TM, NamedParameters>                                  VPM_helper;
  typedef typename VPM_helper::const_type                                                VPM;
  VPM vpmap = VPM_helper::get_const_map(np, tmesh);

  const bool do_limit = !(is_default_parameter<NamedParameters, internal_np::maximum_number_t>::value);
  const unsigned int maximum_number = choose_parameter(get_parameter(np, internal_np::maximum_number), 0);
  if(do_limit && maximum_number == 0)
  {
    return out;
  }
  unsigned int counter = 0;
  const unsigned int seed = choose_parameter(get_parameter(np, internal_np::random_seed), 0);
  CGAL_USE(seed); // used in the random shuffle of the range, which is only done to balance tasks in parallel

  const std::ptrdiff_t cutoff = 2000;

  // make one box per face
  std::vector<Box> boxes;
  boxes.reserve(std::distance(std::begin(face_range), std::end(face_range)));

  // This loop is very cheap, so there is hardly anything to gain from parallelizing it
  for(face_descriptor f : face_range)
  {
    std::array<vertex_descriptor, 3> vh;
    Wrapper::get_face_vertices(f, vh, tmesh);

    typename boost::property_traits<VPM>::reference
      p = get(vpmap, vh[0]),
      q = get(vpmap, vh[1]),
      r = get(vpmap, vh[2]);

    // tiny fixme: if f is degenerate, we might still have a real intersection between f
    // and another face f', but right now we are not creating a box for f and thus not returning those
    if(collinear(p, q, r))
    {
      if(throw_on_SI)
        throw CGAL::internal::Throw_at_output_exception();
      else
      {
        *out++= std::make_pair(f, f);
        ++counter;
        if(do_limit && counter == maximum_number)
        {
          return out;
        }
      }
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
  static_assert (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    // We are going to split the range into a number of smaller ranges. To handle
    // smaller trees of roughly the same size, we first apply a random shuffle to the range
    CGAL::Random rng(seed);
    CGAL::cpp98::random_shuffle(box_ptr.begin(), box_ptr.end(), rng);

    // Write in a concurrent vector all pairs that intersect
    typedef tbb::concurrent_vector<std::pair<face_descriptor, face_descriptor> >         Face_pairs;
    typedef std::back_insert_iterator<Face_pairs>                                        Face_pairs_back_inserter;
    typedef internal::Strict_intersect_faces<Box, TM, VPM, GT, Face_pairs_back_inserter> Intersecting_faces_filter;
    //for maximum_number
    typedef internal::Throw_at_count_reached_functor<Face_pairs_back_inserter>           Throw_functor;
    typedef boost::function_output_iterator<Throw_functor>                               Throwing_after_count_output_iterator;
    typedef internal::Strict_intersect_faces<Box, TM, VPM,
        GT,Throwing_after_count_output_iterator>                                         Filtered_intersecting_faces_filter;

    Face_pairs face_pairs;
    if(throw_on_SI)
      CGAL::box_self_intersection_d<ConcurrencyTag>(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
    else if(do_limit)
    {
      try
      {
        std::atomic<unsigned int> atomic_counter(counter);
        Throw_functor throwing_count_functor(atomic_counter, maximum_number, std::back_inserter(face_pairs));
        Throwing_after_count_output_iterator count_filter(throwing_count_functor);
         Filtered_intersecting_faces_filter limited_callback(tmesh, vpmap, gt, count_filter);
        CGAL::box_self_intersection_d<ConcurrencyTag>(box_ptr.begin(), box_ptr.end(), limited_callback, cutoff);
      }
      catch(const CGAL::internal::Throw_at_output_exception&)
      {
        // Sequentially write into the output iterator
        std::copy(face_pairs.begin(), face_pairs.end(), out);
        return out;
      }
    }
    else
    {
      Intersecting_faces_filter callback(tmesh, vpmap, gt, std::back_inserter(face_pairs));
      CGAL::box_self_intersection_d<ConcurrencyTag>(box_ptr.begin(), box_ptr.end(), callback, cutoff);
    }

    // Sequentially write into the output iterator
    for(std::size_t i=0; i<face_pairs.size(); ++i)
      *out ++= face_pairs[i];

    return out;
  }
#endif

  // Sequential version of the code
  // Compute self-intersections filtered out by boxes
  typedef internal::Strict_intersect_faces<Box, TM, VPM, GT, FacePairOutputIterator> Intersecting_faces_filter;
  Intersecting_faces_filter intersect_faces(tmesh, vpmap, gt, out);

  if(throw_on_SI)
    CGAL::box_self_intersection_d<CGAL::Sequential_tag>(box_ptr.begin(), box_ptr.end(), throwing_filter, cutoff);
  else if(do_limit)
  {
    typedef std::function<void(const std::pair<face_descriptor, face_descriptor>&) > Count_and_throw_filter;
    std::size_t nbi=0;
    Count_and_throw_filter max_inter_counter = [&nbi, maximum_number, &out](const std::pair<face_descriptor, face_descriptor>& f_pair)
    {
      *out++=f_pair;
      if (++nbi == maximum_number)
        throw CGAL::internal::Throw_at_output_exception();
    };
    typedef internal::Strict_intersect_faces<Box, TM, VPM, GT, boost::function_output_iterator<Count_and_throw_filter > > Intersecting_faces_limited_filter;
    Intersecting_faces_limited_filter limited_intersect_faces(tmesh, vpmap, gt,
                                                      boost::make_function_output_iterator(max_inter_counter));
    try
    {
      CGAL::box_self_intersection_d<CGAL::Sequential_tag>(box_ptr.begin(), box_ptr.end(), limited_intersect_faces, cutoff);
    }
    catch (const CGAL::internal::Throw_at_output_exception&)
    {
      return out;
    }
  }
  else
    CGAL::box_self_intersection_d<CGAL::Sequential_tag>(box_ptr.begin(), box_ptr.end(), intersect_faces, cutoff);

  return intersect_faces.m_iterator;
}

} // namespace internal

/*!
 * \ingroup PMP_intersection_grp
 *
 * collects intersections between a subset of faces of a triangulated surface mesh.
 * Two faces are said to intersect if the corresponding triangles intersect
 * and the intersection is neither an edge nor a vertex incident to both faces.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam FaceRange a model of `ConstRange` with value type `boost::graph_traits<TriangleMesh>::%face_descriptor`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePairOutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`.
 *    It does not need to be thread-safe.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param face_range the range of faces to check for self-intersection.
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_number}
 *     \cgalParamDescription{the maximum number of self intersections that will be detected and returned by the function.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{No limit.}
 *     \cgalParamExtra{In parallel mode, the number of returned self-intersections is at least `maximum_number`
 *     (and not exactly that number) as no strong synchronization is put on threads for performance reasons.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `out`
 *
 * @sa `does_self_intersect()`
 */
template < class ConcurrencyTag = Sequential_tag,
           class TriangleMesh,
           class FaceRange,
           class FacePairOutputIterator,
           class NamedParameters = parameters::Default_named_parameters>
FacePairOutputIterator
self_intersections(const FaceRange& face_range,
                   const TriangleMesh& tmesh,
                         FacePairOutputIterator out,
                   const NamedParameters& np = parameters::default_values())
{
  return internal::self_intersections_impl<ConcurrencyTag>(face_range, tmesh, out, false /*don't throw*/, np);
}

/**
 * \ingroup PMP_intersection_grp
 *
 * collects intersections between all the faces of a triangulated surface mesh.
 * Two faces are said to intersect if the corresponding triangles intersect
 * and the intersection is neither an edge nor a vertex incident to both faces.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                         Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam FacePairOutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<boost::graph_traits<TriangleMesh>::%face_descriptor, boost::graph_traits<TriangleMesh>::%face_descriptor>`.
 *    It does not need to be thread-safe.
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tmesh the triangulated surface mesh to be checked
 * @param out output iterator to be filled with all pairs of non-adjacent faces that intersect.
 *            In case `tmesh` contains some degenerate faces, for each degenerate face `f` a pair `(f,f)`
 *            will be put in `out` before any other self intersection between non-degenerate faces. <br>
 *            Note that these are the only pairs where degenerate faces will be reported.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{maximum_number}
 *     \cgalParamDescription{the maximum number of self intersections that will be detected and returned by the function.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{No limit.}
 *     \cgalParamExtra{In parallel mode, the number of returned self-intersections is at least `maximum_number`
 *     (and not exactly that number) as no strong synchronization is put on threads for performance reasons.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `out`
 *
 * @sa `does_self_intersect()`
 */
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh,
          class FacePairOutputIterator,
          class CGAL_NP_TEMPLATE_PARAMETERS>
FacePairOutputIterator
self_intersections(const TriangleMesh& tmesh,
                         FacePairOutputIterator out,
                   const CGAL_NP_CLASS& np = parameters::default_values())
{
  return self_intersections<ConcurrencyTag>(faces(tmesh), tmesh, out, np);
}

/**
 * \ingroup PMP_intersection_grp
 *
 * \brief tests if a set of faces of a triangulated surface mesh self-intersects.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam FaceRange a range of `face_descriptor`
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param face_range the set of faces to test for self-intersection
 * @param tmesh the triangulated surface mesh to be tested
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `true` if the faces in `face_range` self-intersect
 *
 * @sa `self_intersections()`
 */
template <class ConcurrencyTag = Sequential_tag,
          class FaceRange,
          class TriangleMesh,
          class NamedParameters = parameters::Default_named_parameters>
bool does_self_intersect(const FaceRange& face_range,
                         const TriangleMesh& tmesh,
                         const NamedParameters& np = parameters::default_values())
{
  try
  {
    CGAL::Emptyset_iterator unused_out;
    internal::self_intersections_impl<ConcurrencyTag>(face_range, tmesh, unused_out, true /*throw*/, np);
  }
  catch (const CGAL::internal::Throw_at_output_exception&)
  {
    return true;
  }
  #if defined(CGAL_LINKED_WITH_TBB) && TBB_USE_CAPTURED_EXCEPTION
  catch (const tbb::captured_exception& e)
  {
    const char* ti1 = e.name();
    const char* ti2 = typeid(const CGAL::internal::Throw_at_output_exception&).name();
    const std::string tn1(ti1);
    const std::string tn2(ti2);
    if (tn1 == tn2) return true;
    else throw;
  }
  #endif
  return false;
}

/**
 * \ingroup PMP_intersection_grp
 *
 * \brief tests if a triangulated surface mesh self-intersects.
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @pre `CGAL::is_triangle_mesh(tmesh)`
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam TriangleMesh a model of `FaceListGraph`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param tmesh the triangulated surface mesh to be tested
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
 *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor` as key type and `%Point_3` as value type}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
 *                     should be available for the vertices of `tmesh`.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `true` if `tmesh` self-intersects
 *
 * @sa `self_intersections()`
 */
template <class ConcurrencyTag = Sequential_tag,
          class TriangleMesh,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool does_self_intersect(const TriangleMesh& tmesh,
                         const CGAL_NP_CLASS& np = parameters::default_values())
{
  return does_self_intersect<ConcurrencyTag>(faces(tmesh), tmesh, np);
}


#ifndef DOXYGEN_RUNNING
template <class PointRange, class VPM>
struct Property_map_for_soup
{
  typedef std::size_t key_type;
  typedef typename boost::property_traits<VPM>::value_type value_type;
  //typedef typename boost::property_traits<VPM>::category category;
  typedef boost::readable_property_map_tag category;
  typedef typename boost::property_traits<VPM>::reference reference;

  const PointRange& points;
  VPM vpm;

  Property_map_for_soup(const PointRange& points, VPM vpm)
    : points(points)
    , vpm(vpm)
  {}

  inline friend
  reference get(const Property_map_for_soup<PointRange, VPM>& map, key_type k)
  {
    return get(map.vpm, map.points[k]);
  }
};
#endif

/**
 * \ingroup PMP_intersection_grp
 *
 * collects intersections between all the triangles in a triangle soup.
 *
 * Two triangles of the soup are said to intersect if the corresponding geometric triangles intersect
 * and the intersection is neither an edge nor a vertex of both triangles
 * (with the same point ids, ignoring the orientation for an edge).
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam PointRange a model of the concept `RandomAccessContainer`
 *         whose value type is the point type
 * @tparam TriangleRange a model of the concept `RandomAccessContainer` whose
 *         value type is a model of the concept `RandomAccessContainer` whose value type is `std::size_t`
 * @tparam TriangleIdPairOutputIterator a model of `OutputIterator` holding objects of type
 *   `std::pair<std::size_t,std::size_t>`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param points points of the soup of triangles
 * @param triangles each element in the range describes a triangle using the indices of the points in `points`
 * @param out output iterator to be filled with all pairs of ids of triangles intersecting (the id of a triangle is its position in `triangles`).
 *            In case the triangle soup contains some degenerate triangles, for each degenerate triangle `t` with id `i` a pair `(i,i)`
 *            will be put in `out` before any other self intersection between non-degenerate triangles.<br>
 *            Note that these are the only pairs where degenerate triangles will be reported.
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{maximum_number}
 *     \cgalParamDescription{the maximum number of self intersections that will be detected and returned by the function.}
 *     \cgalParamType{unsigned int}
 *     \cgalParamDefault{No limit.}
 *     \cgalParamExtra{In parallel mode, the number of returned self-intersections is at least `maximum_number`
 *     (and not exactly that number) as no strong synchronization is put on threads for performance reasons.}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is a point type from a \cgal `Kernel`.}
 *     \cgalParamDefault{`CGAL::Identity_property_map`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the point type of the point map.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `out`
 *
 * @sa `does_triangle_soup_self_intersect()`
 * @sa `self_intersections()`
 * @sa `does_self_intersect()`
 */
template <class ConcurrencyTag = Sequential_tag,
          class PointRange,
          class TriangleRange,
          class TriangleIdPairOutputIterator,
          class CGAL_NP_TEMPLATE_PARAMETERS>
TriangleIdPairOutputIterator
triangle_soup_self_intersections(const PointRange& points,
                                 const TriangleRange& triangles,
                                 TriangleIdPairOutputIterator out,
                                 const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::const_type Point_map_base;
  Point_map_base pm_base = choose_parameter<Point_map_base>(get_parameter(np, internal_np::point_map));
  typedef Property_map_for_soup<PointRange, Point_map_base> Point_map;
  typedef typename GetPolygonSoupGeomTraits<PointRange, CGAL_NP_CLASS>::type GT;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  const bool do_limit = !(is_default_parameter<CGAL_NP_CLASS, internal_np::maximum_number_t>::value);
  if (do_limit)
  {
    return self_intersections<ConcurrencyTag>(boost::irange<std::size_t>(0, triangles.size()),
                                              std::make_pair(std::cref(points), std::cref(triangles)),
                                              out,
                                              parameters::vertex_point_map(Point_map(points,pm_base)).
                                              geom_traits(gt).
                                              maximum_number(choose_parameter(get_parameter(np, internal_np::maximum_number), 0)));
  }

  return self_intersections<ConcurrencyTag>(boost::irange<std::size_t>(0, triangles.size()),
                                            std::make_pair(std::cref(points), std::cref(triangles)),
                                            out,
                                            parameters::vertex_point_map(Point_map(points,pm_base)).
                                            geom_traits(gt));
}

/**
 * \ingroup PMP_intersection_grp
 *
 * \brief tests if a triangle soup self-intersects.
 *
 * A triangle soup self-intersects if at least two triangles of the soup intersect.
 * Two triangles of the soup are said to intersect if the corresponding geometric triangles intersect
 * and the intersection is neither an edge nor a vertex of both triangles
 * (with the same point ids, ignoring the orientation for an edge).
 *
 * This function depends on the package \ref PkgBoxIntersectionD.
 *
 * @tparam ConcurrencyTag enables sequential versus parallel algorithm.
 *                        Possible values are `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * @tparam PointRange a model of the concept `RandomAccessContainer`
 *         whose value type is the point type
 * @tparam TriangleRange a model of the concept `RandomAccessContainer` whose
 *         value type is a model of the concept `RandomAccessContainer` whose value type is `std::size_t`
 * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * @param points points of the soup of triangles
 * @param triangles each element in the range describes a triangle using the indices of the points in `points`
 * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
 *
 * \cgalNamedParamsBegin
 *   \cgalParamNBegin{point_map}
 *     \cgalParamDescription{a property map associating points to the elements of the range `points`}
 *     \cgalParamType{a model of `ReadablePropertyMap` whose value type is a point type from a \cgal `Kernel`.}
 *     \cgalParamDefault{`CGAL::Identity_property_map`}
 *   \cgalParamNEnd
 *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `PMPSelfIntersectionTraits`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the point type of the point map.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * @return `true` if the triangle soup self-intersects, and `false` otherwise.
 *
 * @sa `triangle_soup_self_intersections()`
 * @sa `self_intersections()`
 * @sa `does_self_intersect()`
 */
template <class ConcurrencyTag = Sequential_tag,
          class PointRange,
          class TriangleRange,
          class CGAL_NP_TEMPLATE_PARAMETERS>
bool does_triangle_soup_self_intersect(const PointRange& points,
                                       const TriangleRange& triangles,
                                       const CGAL_NP_CLASS& np = parameters::default_values())
{
  try
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    CGAL::Emptyset_iterator unused_out;
    typedef typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::const_type Point_map_base;
    Point_map_base pm_base = choose_parameter<Point_map_base>(get_parameter(np, internal_np::point_map));
    typedef Property_map_for_soup<PointRange, Point_map_base> Point_map;
    typedef typename GetPolygonSoupGeomTraits<PointRange, CGAL_NP_CLASS>::type GT;
    GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

    internal::self_intersections_impl<ConcurrencyTag>(boost::irange<std::size_t>(0, triangles.size()),
                                                      std::make_pair(std::cref(points), std::cref(triangles)),
                                                      unused_out, true /*throw*/,
                                                      parameters::vertex_point_map(Point_map(points,pm_base))
                                                                 .geom_traits(gt));
  }
  catch (const CGAL::internal::Throw_at_output_exception&)
  {
    return true;
  }
  #if defined(CGAL_LINKED_WITH_TBB) && TBB_USE_CAPTURED_EXCEPTION
  catch (const tbb::captured_exception& e)
  {
    const char* ti1 = e.name();
    const char* ti2 = typeid(const CGAL::internal::Throw_at_output_exception&).name();
    const std::string tn1(ti1);
    const std::string tn2(ti2);
    if (tn1 == tn2) return true;
    else throw;
  }
  #endif
  return false;
}

}// namespace Polygon_mesh_processing
}// namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_SELF_INTERSECTIONS
