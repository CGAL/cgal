// Copyright (c) 2015, 2018 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     :  Konstantinos Katrioplas,
//                  Mael Rouxel-Labbé

#ifndef CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H
#define CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H

#include <CGAL/license/Polygon_mesh_processing/repair.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#include <CGAL/array.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_kernel_selector.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/range/has_range_iterator.hpp>
#include <boost/graph/graph_traits.hpp>

#include <array>
#include <limits>
#include <map>
#include <utility>
#include <vector>

namespace CGAL {

namespace Polygon_mesh_processing {

/// \ingroup PMP_predicates_grp
///
/// checks whether an edge is degenerate.
/// An edge is considered degenerate if the geometric positions of its two extremities are identical.
///
/// @tparam PolygonMesh a model of `HalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param e an edge of `pm`
/// @param pm polygon mesh containing `e`
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `pm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, pm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functor `Equal_3` to check whether two points are identical.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return `true` if the edge `e` is degenerate, `false` otherwise.
///
/// \sa `degenerate_edges()`
template <typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
bool is_degenerate_edge(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
                        const PolygonMesh& pm,
                        const NamedParameters& np = parameters::default_values())
{
  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_const_property_map(vertex_point, pm));

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Traits;
  Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

  return traits.equal_3_object()(get(vpmap, source(e, pm)), get(vpmap, target(e, pm)));
}

/// \ingroup PMP_predicates_grp
///
/// collects the degenerate edges within a given range of edges.
///
/// @tparam EdgeRange a model of `Range` with value type `boost::graph_traits<TriangleMesh>::%edge_descriptor`
/// @tparam TriangleMesh a model of `HalfedgeGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param edges a subset of edges of `tm`
/// @param tm a triangle mesh
/// @param out an output iterator in which the degenerate edges are written
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `TriangleMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functor `Equal_3` to check whether two points are identical.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \sa `is_degenerate_edge()`
// \sa `remove_degenerate_edges()`
template <typename EdgeRange, typename TriangleMesh, typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator degenerate_edges(const EdgeRange& edges,
                                const TriangleMesh& tm,
                                OutputIterator out,
                                const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::edge_descriptor edge_descriptor;

  for(edge_descriptor ed : edges)
    if(is_degenerate_edge(ed, tm, np))
      *out++ = ed;

  return out;
}

/// \ingroup PMP_predicates_grp
///
/// calls the function `degenerate_edges()` with the range: `edges(tm)`.
///
/// See the other overload for the comprehensive description of the parameters.
template <typename TriangleMesh, typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator degenerate_edges(const TriangleMesh& tm,
                                OutputIterator out,
                                const CGAL_NP_CLASS& np = parameters::default_values()
                               )
{
  return degenerate_edges(edges(tm), tm, out, np);
}

/// \ingroup PMP_predicates_grp
///
/// checks whether a triangle face is degenerate.
/// A triangle face is considered degenerate if the geometric positions of its vertices are collinear.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm a triangle mesh containing `f`
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functor `Collinear_3` to check whether three points are aligned.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return `true` if the face `f` is degenerate, `false` otherwise.
///
/// \sa `degenerate_faces()`
template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
bool is_degenerate_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                                 const TriangleMesh& tm,
                                 const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(CGAL::is_triangle(halfedge(f, tm), tm));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type Traits;
  Traits traits = choose_parameter<Traits>(get_parameter(np, internal_np::geom_traits));

  typename boost::graph_traits<TriangleMesh>::halfedge_descriptor h = halfedge(f, tm);

  return traits.collinear_3_object()(get(vpmap, source(h, tm)),
                                     get(vpmap, target(h, tm)),
                                     get(vpmap, target(next(h, tm), tm)));
}

/// \ingroup PMP_predicates_grp
///
/// collects the degenerate faces within a given range of faces.
///
/// @tparam FaceRange a model of `Range` with value type `boost::graph_traits<TriangleMesh>::%face_descriptor`
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param faces a subset of faces of `tm`
/// @param tm a triangle mesh
/// @param out an output iterator in which the degenerate faces are put
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
///                     must be available in `TriangleMesh`.}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    and the nested functor `Collinear_3` to check whether three points are collinear.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \sa `is_degenerate_triangle_face()`
// `\sa remove_degenerate_faces()`
template <typename FaceRange, typename TriangleMesh, typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator degenerate_faces(const FaceRange& faces,
                                const TriangleMesh& tm,
                                OutputIterator out,
                                const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  for(face_descriptor fd : faces)
  {
    if(is_degenerate_triangle_face(fd, tm, np))
      *out++ = fd;
  }
  return out;
}

/// \ingroup PMP_predicates_grp
///
/// calls the function `degenerate_faces()` with the range: `faces(tm)`.
///
/// See the other overload for the comprehensive description of the parameters.
template <typename TriangleMesh, typename OutputIterator, typename CGAL_NP_TEMPLATE_PARAMETERS>
OutputIterator degenerate_faces(const TriangleMesh& tm,
                                OutputIterator out,
                                const CGAL_NP_CLASS& np = parameters::default_values()
                                )
{
  return degenerate_faces(faces(tm), tm, out, np);
}

namespace internal {

template <typename K, template <class Kernel> class Pred>
struct Get_filtered_predicate_RT
{
  typedef typename Exact_kernel_selector<K>::Exact_kernel_rt                   Exact_kernel_rt;
  typedef typename Exact_kernel_selector<K>::C2E_rt                            C2E_rt;
  typedef Simple_cartesian<Interval_nt_advanced>                               Approximate_kernel;
  typedef Cartesian_converter<K, Approximate_kernel>                           C2F;

  typedef Filtered_predicate<Pred<Exact_kernel_rt>,
                             Pred<Approximate_kernel>,
                             C2E_rt, C2F> type;
};

// predicates

template <typename K>
struct Is_edge_length_ratio_over_threshold_impl
{
  typedef int result_type;

  /// Computes the ratio between the longest edge's length and the shortest edge's length
  /// and compares with a user-defined bound.
  /// Returns -1 if the ratio is below the bound, and 0, 1, or 2 otherwise, with the value
  /// indicating the shortest halfedge.
  int operator()(const typename K::Point_3& p,
                 const typename K::Point_3& q,
                 const typename K::Point_3& r,
                 const typename K::FT threshold_squared) const
  {
    typedef typename K::FT                                                     FT;

    FT sq_length_0 = K().compute_squared_distance_3_object()(p, q);

    FT min_sq_length = sq_length_0, max_sq_length = sq_length_0;
    int min_id = 0;

    auto get_min_max = [&](const typename K::Point_3& pa, const typename K::Point_3& pb, int id) -> void
    {
      const FT sq_length = K().compute_squared_distance_3_object()(pa, pb);

      if(max_sq_length < sq_length)
        max_sq_length = sq_length;

      if(sq_length < min_sq_length)
      {
        min_sq_length = sq_length;
        min_id = id;
      }
    };

    get_min_max(q, r, 1);
    get_min_max(r, p, 2);

    if(min_sq_length == 0)
      return min_id;

    if(compare(max_sq_length, threshold_squared * min_sq_length) != CGAL::SMALLER)
      return min_id;
    else
      return -1;
  }
};

template<typename K, bool has_filtered_predicates = K::Has_filtered_predicates>
struct Is_edge_length_ratio_over_threshold
  : public Is_edge_length_ratio_over_threshold_impl<K>
{
  using Is_edge_length_ratio_over_threshold_impl<K>::operator();
};

template<typename K>
struct Is_edge_length_ratio_over_threshold<K, true>
  : public Get_filtered_predicate_RT<K, Is_edge_length_ratio_over_threshold_impl>::type
{
  using Get_filtered_predicate_RT<K, Is_edge_length_ratio_over_threshold_impl>::type::operator();
};

} // namespace internal

/// \ingroup PMP_predicates_grp
///
/// checks whether a triangle face is needle.
/// A triangle is said to be a <i>needle</i> if its longest edge is much longer than its shortest edge.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm triangle mesh containing `f`
/// @param threshold a bound on the ratio of the longest edge length and the shortest edge length
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `FT`,
///                    and the nested functor `Compute_squared_distance_3`.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the shortest halfedge if the triangle face is a needle, and a null halfedge otherwise.
///         If the face contains degenerate edges, a halfedge corresponding to one of these edges is returned.
///
/// \sa `is_cap_triangle_face()`
template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_needle_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                        const TriangleMesh& tm,
                        const double threshold,
                        const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(threshold >= 1.);
  CGAL_precondition(f != boost::graph_traits<TriangleMesh>::null_face());
  CGAL_precondition(CGAL::is_triangle(halfedge(f, tm), tm));

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor       halfedge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type           Traits;

  const halfedge_descriptor h = halfedge(f, tm);

  internal::Is_edge_length_ratio_over_threshold<Traits> pred;
  const int res = pred(get(vpmap, source(h, tm)),
                       get(vpmap, target(h, tm)),
                       get(vpmap, target(next(h, tm), tm)),
                       square(threshold));

  if(res == -1)
    return boost::graph_traits<TriangleMesh>::null_halfedge();
  if(res == 0)
    return h;
  if(res == 1)
    return next(h, tm);
  else
    return prev(h, tm);
}

namespace internal {

template <typename K>
struct Is_cap_angle_over_threshold_impl
{
  typedef int result_type;

  /// Computes the ratio between the longest edge's length and the shortest edge's length
  /// and compares with a user-defined bound.
  /// Returns -1 if the ratio is below the bound, and 0, 1, or 2 otherwise, with the value
  /// indicating the shortest halfedge.
  int operator()(const typename K::Point_3& p,
                 const typename K::Point_3& q,
                 const typename K::Point_3& r,
                 const typename K::FT threshold_squared) const
  {
    typedef typename K::FT                                                     FT;
    typedef typename K::Vector_3                                               Vector_3;

    std::array<FT, 3> sq_lengths = { K().compute_squared_distance_3_object()(p, q),
                                     K().compute_squared_distance_3_object()(q, r),
                                     K().compute_squared_distance_3_object()(r, p) };

    // If even one edge is degenerate, it cannot be a cap
    if(is_zero(sq_lengths[0]) || is_zero(sq_lengths[1]) || is_zero(sq_lengths[2]))
      return -1;

    auto handle_triplet = [&](const typename K::Point_3& pa,
                              const typename K::Point_3& pb,
                              const typename K::Point_3& pc, int pos) -> bool
    {
      const Vector_3 vc = K().construct_vector_3_object()(pb, pc);
      const Vector_3 va = K().construct_vector_3_object()(pb, pa);
      const FT dot_ca = K().compute_scalar_product_3_object()(vc, va);
      const bool neg_sp = !(is_positive(dot_ca));
      if(!neg_sp)
        return false;

      const FT sq_c = sq_lengths[(pos+1)%3];
      const FT sq_a = sq_lengths[pos];

      return (compare(square(dot_ca), threshold_squared * sq_c * sq_a) != CGAL::SMALLER);
    };

    // halfedge 0 is between p and q, so cap at q => return halfedge 2 (r to p)
    if(handle_triplet(p, q, r, 0))
      return 2;
    if(handle_triplet(q, r, p, 1))
      return 0;
    if(handle_triplet(r, p, q, 2))
      return 1;

    return -1;
  }
};

template<typename K, bool has_filtered_predicates = K::Has_filtered_predicates>
struct Is_cap_angle_over_threshold
  : public Is_cap_angle_over_threshold_impl<K>
{
  using Is_cap_angle_over_threshold_impl<K>::operator();
};

template<typename K>
struct Is_cap_angle_over_threshold<K, true>
  : public Get_filtered_predicate_RT<K, Is_cap_angle_over_threshold_impl>::type
{
  using Get_filtered_predicate_RT<K, Is_cap_angle_over_threshold_impl>::type::operator();
};

} // namespace internal

/// \ingroup PMP_predicates_grp
///
/// checks whether a triangle face is a cap.
/// A triangle is said to be a <i>cap</i> if one of the its angles is close to `180` degrees.
///
/// @tparam TriangleMesh a model of `FaceGraph`
/// @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
///
/// @param f a triangle face of `tm`
/// @param tm triangle mesh containing `f`
/// @param threshold the cosine of a minimum angle such that if `f` has an angle greater than this bound,
///                  it is a cap. The threshold is in range `[-1 0]` and corresponds to an angle
///                  between `90` and `180` degrees.
/// @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
///
/// \cgalNamedParamsBegin
///   \cgalParamNBegin{vertex_point_map}
///     \cgalParamDescription{a property map associating points to the vertices of `tm`}
///     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
///                    as key type and `%Point_3` as value type}
///     \cgalParamDefault{`boost::get(CGAL::vertex_point, tm)`}
///   \cgalParamNEnd
///
///   \cgalParamNBegin{geom_traits}
///     \cgalParamDescription{an instance of a geometric traits class}
///     \cgalParamType{The traits class must provide the nested type `Point_3`,
///                    the nested functors `Compute_squared_distance_3`, `Construct_vector_3`,
///                    and `Compute_scalar_product_3`.}
///     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
///     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
///   \cgalParamNEnd
/// \cgalNamedParamsEnd
///
/// \return the halfedge opposite of the largest angle if the face is a cap, and a null halfedge otherwise.
///
/// \sa `is_needle_triangle_face()`
template <typename TriangleMesh, typename NamedParameters = parameters::Default_named_parameters>
typename boost::graph_traits<TriangleMesh>::halfedge_descriptor
is_cap_triangle_face(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                     const TriangleMesh& tm,
                     const double threshold,
                     const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(f != boost::graph_traits<TriangleMesh>::null_face());
  CGAL_precondition(CGAL::is_triangle(halfedge(f, tm), tm));
  CGAL_precondition(threshold >= -1.);
  CGAL_precondition(threshold <= 0.);

  using parameters::get_parameter;
  using parameters::choose_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor   halfedge_descriptor;

  typedef typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type VertexPointMap;
  VertexPointMap vpmap = choose_parameter(get_parameter(np, internal_np::vertex_point),
                                          get_const_property_map(vertex_point, tm));

  typedef typename GetGeomTraits<TriangleMesh, NamedParameters>::type       Traits;

  const halfedge_descriptor h = halfedge(f, tm);

  internal::Is_cap_angle_over_threshold<Traits> pred;
  const int res = pred(get(vpmap, source(h, tm)),
                       get(vpmap, target(h, tm)),
                       get(vpmap, target(next(h, tm), tm)),
                       square(threshold));

  if(res == -1)
    return boost::graph_traits<TriangleMesh>::null_halfedge();
  if(res == 0)
    return h;
  if(res == 1)
    return next(h, tm);
  else
    return prev(h, tm);
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_SHAPE_PREDICATES_H
