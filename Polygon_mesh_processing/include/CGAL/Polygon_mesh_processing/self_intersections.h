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

#include <CGAL/Polygon_mesh_processing/internal/soup_and_face_graph_tree_helper.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_indexed_triangle_primitive_3.h>
#include <CGAL/AABB_trees/intersection.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <class TM, class GT, class VPM>
struct Triangle_mesh_and_triangle_soup_wrapper
{
  using face_descriptor     = typename boost::graph_traits<TM>::face_descriptor;
  using vertex_descriptor   = typename boost::graph_traits<TM>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TM>::halfedge_descriptor; // private

  using Tree_helper = AABB_tree_graph_helper<TM, GT, VPM>;
  using Tree        = typename Tree_helper::Tree;

  template<class ConcurrencyTag = Sequential_tag, class FaceRange>
  static void build_tree(const FaceRange& faces, Tree& tree, const TM &tm, VPM vpm)
  {
    Tree_helper().template build<ConcurrencyTag>(faces, tree, tm, vpm);
  }

  template<class ConcurrencyTag = Sequential_tag>
  static void build_tree(Tree& tree, const TM &tm, VPM vpm)
  {
    Tree_helper().template build<ConcurrencyTag>(tree, tm, vpm);
  }

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

  static constexpr bool allow_identical_face()
  {
    return false;
  }
};

template <class PointRange, class TriangleRange, class GT, class VPM>
struct Triangle_mesh_and_triangle_soup_wrapper< std::pair<const PointRange&, const TriangleRange&>, GT, VPM>
{
  using face_descriptor  = std::size_t;
  using vertex_descriptor = std::size_t;

  using Soup = std::pair<const PointRange&, const TriangleRange& >;

  using Tree_helper = AABB_tree_soup_helper<PointRange, TriangleRange, GT, VPM>;
  using Tree        = typename Tree_helper::Tree;

  template<class ConcurrencyTag=Sequential_tag, class FaceRange>
  static void build_tree(const FaceRange &faces, Tree& tree, const Soup& soup, VPM vpm){
    Tree_helper().template build<ConcurrencyTag>(faces, tree, soup.first, soup.second, vpm);
  }

  template<class ConcurrencyTag=Sequential_tag>
  static void build_tree(Tree& tree, const Soup& soup, VPM vpm){
    Tree_helper().template build<ConcurrencyTag>( tree, soup.first, soup.second, vpm);
  }

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

  static constexpr bool allow_identical_face()
  {
    return true;
  }
};

template<typename OutputIterator>
struct Throw_at_count_reached_output_iterator
{
  using Self = Throw_at_count_reached_output_iterator<OutputIterator>;
  std::atomic<unsigned int> &counter;
  const unsigned int &maxval;
  OutputIterator out;

  using iterator_category = std::output_iterator_tag;
  Throw_at_count_reached_output_iterator(std::atomic<unsigned int> &counter,
                                         const unsigned int &maxval,
                                         OutputIterator out)
    : counter(counter), maxval(maxval), out(out)
  {}

  template<class T>
  Self& operator=(const T& t)
  {
    *out++ = t;
    ++counter;
    if(counter >= maxval)
      throw CGAL::internal::Throw_at_output_exception();
    return *this;
  }

  Self& operator*(){ return *this; }
  Self& operator++(){ return *this; }
  Self& operator++(int){ return *this; }
};

// Checks for 'real' intersections, i.e. not simply a shared vertex or edge
template <class GT, class TM, class VPM>
bool do_faces_intersect(typename Triangle_mesh_and_triangle_soup_wrapper<TM, GT, VPM>::face_descriptor fh,
                        typename Triangle_mesh_and_triangle_soup_wrapper<TM, GT, VPM>::face_descriptor fg,
                        const TM& tmesh,
                        const VPM vpmap,
                        const typename GT::Construct_segment_3& construct_segment,
                        const typename GT::Construct_triangle_3& construct_triangle,
                        const typename GT::Do_intersect_3& do_intersect)
{
  using Wrapper = Triangle_mesh_and_triangle_soup_wrapper<TM, GT, VPM>;
  using vertex_descriptor = typename Wrapper::vertex_descriptor;

  using Segment  = typename GT::Segment_3;
  using Triangle = typename GT::Triangle_3;

  std::array<vertex_descriptor, 3> hv, gv;
  Wrapper::get_face_vertices(fh, hv, tmesh);
  Wrapper::get_face_vertices(fg, gv, tmesh);

  // check for shared edge
  std::array<vertex_descriptor, 4> verts;
  if (Wrapper::faces_have_a_shared_edge(fh, fg, verts, tmesh))
  {
    if (Wrapper::allow_identical_face() && verts[2]==verts[3]) return false; // only for a soup of triangles

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

template<typename AABBTraits, typename GT, typename OutputIterator, typename TriangleMesh, typename VertexPointMap>
class Listing_distinct_intersecting_faces_traits
{
  using FT = typename AABBTraits::FT;
  using Point = typename AABBTraits::Point;
  using Primitive = typename AABBTraits::Primitive;
  using Bounding_box = typename AABBTraits::Bounding_box;
  using Primitive_id = typename AABBTraits::Primitive::Id;
  using Point_and_primitive_id = typename AABBTraits::Point_and_primitive_id;
  using Object_and_primitive_id = typename AABBTraits::Object_and_primitive_id;
  using Node = ::CGAL::AABB_node<AABBTraits>;

public:
  Listing_distinct_intersecting_faces_traits(OutputIterator out_it, const GT& gt, const AABBTraits& traits, const TriangleMesh &tm, const VertexPointMap &vpm)
    : m_out_it(out_it), m_gt(gt), m_traits(traits), m_tm(tm), m_vpm(vpm) {}

  constexpr bool go_further() const { return true; }
  void intersection(const Primitive& query, const Primitive& primitive)
  {
    if(do_faces_intersect<GT>(query.id(), primitive.id(), m_tm, m_vpm,
                              m_gt.construct_segment_3_object(),
                              m_gt.construct_triangle_3_object(),
                              m_gt.do_intersect_3_object()))
    {
      *m_out_it++ = std::make_pair(query.id(), primitive.id());
    }
  }

  bool do_intersect(const Primitive& query, const Node& node) const
  {
    return m_traits.do_intersect_object()(CGAL::internal::Primitive_helper<AABBTraits>::get_datum(query, m_traits).bbox(), node.bbox());
  }

private:
  OutputIterator m_out_it;
  const GT& m_gt;
  const AABBTraits& m_traits;
  const TriangleMesh& m_tm;
  const VertexPointMap& m_vpm;
};

template<typename AABBTraits, typename GT, typename OutputIterator, typename TriangleMesh, typename VertexPointMap>
class Listing_intersecting_faces_two_trees_traits
{
  using Primitive = typename AABBTraits::Primitive;
  using Node = ::CGAL::AABB_node<AABBTraits>;

public:
  Listing_intersecting_faces_two_trees_traits(const AABBTraits& traits,
                                              const GT& gt,
                                              OutputIterator out_,
                                              const TriangleMesh& tm,
                                              const VertexPointMap& vpm)
    : m_traits(traits), m_gt(gt), out(out_), m_tm(tm), m_vpm(vpm) {}

    constexpr bool go_further() const { return true; }

  // node coming from the same tree, always alternate is a good choice
  template<class Node_A, class Node_B>
  bool prefer_A_for_next_step(const Node_A&, const Node_B&, const std::size_t&, const std::size_t&) const {
    return false;
  }

  void intersection(const Primitive& primitive1, const Primitive& primitive2)
  {
    if(do_faces_intersect<GT>(primitive1.id(), primitive2.id(), m_tm, m_vpm,
                              m_gt.construct_segment_3_object(),
                              m_gt.construct_triangle_3_object(),
                              m_gt.do_intersect_3_object()))
    {
      *out++ = std::make_pair(primitive1.id(), primitive2.id());
    }
  }

  void intersection(const Primitive& primitive1, const Node& node2, std::size_t nb_primitives_2)
  {
    Listing_distinct_intersecting_faces_traits<AABBTraits, GT, OutputIterator, TriangleMesh, VertexPointMap> traits(out, m_gt, m_traits, m_tm, m_vpm);
    node2.traversal( primitive1, traits, nb_primitives_2);
  }

  void intersection(const Node& node1, std::size_t nb_primitives_1, const Primitive& primitive2)
  {
    Listing_distinct_intersecting_faces_traits<AABBTraits, GT, OutputIterator, TriangleMesh, VertexPointMap> traits(out, m_gt, m_traits, m_tm, m_vpm);
    node1.traversal( primitive2, traits, nb_primitives_1);
  }

  bool do_intersect(const Node& node1, const Node& node2) const
  {
    return do_overlap(node1.bbox(), node2.bbox());
  }

private:
  const AABBTraits& m_traits;
  const GT& m_gt;
  OutputIterator out;
  const TriangleMesh& m_tm;
  const VertexPointMap& m_vpm;
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
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;
  using CGAL::parameters::is_default_parameter;

  using TM = TriangleMesh;
  using GT = typename GetGeomTraits<TM, NamedParameters>::type;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using VPM_helper = GetVertexPointMap<TM, NamedParameters>;
  using VPM = typename VPM_helper::const_type;
  VPM vpmap = VPM_helper::get_const_map(np, tmesh);

  const bool do_limit = !(is_default_parameter<NamedParameters, internal_np::maximum_number_t>::value);
  const unsigned int maximum_number = choose_parameter(get_parameter(np, internal_np::maximum_number), 0);
  if(do_limit && maximum_number == 0)
  {
    return out;
  }
  unsigned int counter = 0;

  using Wrapper = Triangle_mesh_and_triangle_soup_wrapper<TM, GT, VPM>;
  using face_descriptor = typename Wrapper::face_descriptor;
  using vertex_descriptor = typename Wrapper::vertex_descriptor;

  using AABB_tree = typename Wrapper::Tree;
  using AABB_traits = typename AABB_tree::AABB_traits;

  CGAL_precondition(Wrapper::is_pure_triangle(tmesh));

  // This loop is very cheap, so there is hardly anything to gain from parallelizing it
  std::vector<face_descriptor> faces_not_degenerated;
  faces_not_degenerated.reserve(face_range.size());
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
    if(collinear(p, q, r)){
      if(throw_on_SI)
        throw CGAL::internal::Throw_at_output_exception();
      else
      {
        *out++= std::make_pair(f, f);
        ++counter;
        if(do_limit && counter == maximum_number)
          return out;
      }
    } else {
      faces_not_degenerated.push_back(f);
    }
  }

  // In case we are throwing, like in `does_self_intersect()`, we keep the geometric test to throw ASAP.
  // This is obviously not optimal if there are no or few self-intersections: it would be a greater speed-up
  // to do the same as for `self_intersections()`. However, doing like `self_intersections()` would
  // be a major slow-down over sequential code if there are a lot of self-intersections...
  using Throwing_output_iterator = boost::function_output_iterator<CGAL::internal::Throw_at_output>;
  Throwing_output_iterator throwing_filter;

  AABB_tree tree;
  Wrapper::template build_tree<ConcurrencyTag>(faces_not_degenerated, tree, tmesh, vpmap);

  const auto all_pairs_of_intersecting_faces = [&](const AABB_tree& tree, auto out){
    Listing_intersecting_faces_two_trees_traits<AABB_traits, GT, decltype(out), TM, VPM>
      traversal_traits(tree.traits(), gt, out, tmesh, vpmap);
    CGAL::internal::AABB_tree::one_tree_traversal<ConcurrencyTag>(tree, traversal_traits);
  };
#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  if(std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    // Write in a concurrent vector all pairs that intersect
    using Face_pairs = tbb::concurrent_vector<std::pair<face_descriptor, face_descriptor> >;
    using Face_pairs_inserter = std::back_insert_iterator<Face_pairs>;
    //for maximum_number
    using Throw_iterator = internal::Throw_at_count_reached_output_iterator<Face_pairs_inserter>;

    Face_pairs face_pairs;
    if(throw_on_SI)
      all_pairs_of_intersecting_faces(tree, throwing_filter);
    else if(do_limit)
    {
      try
      {
        std::atomic<unsigned int> atomic_counter(counter);
        Throw_iterator throwing_count(atomic_counter, maximum_number, std::back_inserter(face_pairs));
        Throw_at_count_reached_output_iterator count_filter(throwing_count);
        all_pairs_of_intersecting_faces(tree, count_filter);
      }
      catch(const CGAL::internal::Throw_at_output_exception&)
      {
        // Sequentially write into the output iterator
        for(std::size_t i=0; i<face_pairs.size(); ++i)
          *out ++= face_pairs[i];
      }
    }
    else
    {
      all_pairs_of_intersecting_faces(tree, std::back_inserter(face_pairs));
    }

    // Sequentially write into the output iterator
    for(std::size_t i=0; i<face_pairs.size(); ++i)
      *out ++= face_pairs[i];
    return out;
  }
  else
#endif

  // Sequential version of the code
  if(throw_on_SI)
    all_pairs_of_intersecting_faces(tree, throwing_filter);
  else if(do_limit)
  {
    using Count_and_throw_filter = std::function<void(const std::pair<face_descriptor, face_descriptor>&) >;
    std::size_t nbi=0;
    Count_and_throw_filter max_inter_counter = [&nbi, maximum_number, &out](const std::pair<face_descriptor, face_descriptor>& f_pair)
    {
      *out++=f_pair;
      if (++nbi == maximum_number)
        throw CGAL::internal::Throw_at_output_exception();
    };

    try
    {
      all_pairs_of_intersecting_faces(tree, boost::make_function_output_iterator(max_inter_counter));
    }
    catch (const CGAL::internal::Throw_at_output_exception&)
    {
      return out;
    }
  }
  else
    all_pairs_of_intersecting_faces(tree, out);
  return out;
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
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tmesh)` \endlink
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
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tmesh)` \endlink
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
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tm)` \endlink
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
 * @pre \link CGAL::is_triangle_mesh `CGAL::is_triangle_mesh(tmesh)` \endlink
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

  using Point_map_base = typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::const_type;
  Point_map_base pm_base = choose_parameter<Point_map_base>(get_parameter(np, internal_np::point_map));
  using Point_map = internal::Property_map_for_soup<PointRange, Point_map_base>;
  using GT = typename GetPolygonSoupGeomTraits<PointRange, CGAL_NP_CLASS>::type;
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
    using Point_map_base = typename CGAL::GetPointMap<PointRange, CGAL_NP_CLASS>::const_type;
    Point_map_base pm_base = choose_parameter<Point_map_base>(get_parameter(np, internal_np::point_map));
    using Point_map = internal::Property_map_for_soup<PointRange, Point_map_base>;
    using GT = typename GetPolygonSoupGeomTraits<PointRange, CGAL_NP_CLASS>::type;
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
