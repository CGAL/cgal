// Copyright (c) 2015 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Andreas Fabri

#ifndef CGAL_POLYGON_MESH_PROCESSING_MEASURE_H
#define CGAL_POLYGON_MESH_PROCESSING_MEASURE_H

#include <CGAL/license/Polygon_mesh_processing/measure.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/for_each.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Lazy.h> // needed for CGAL::exact(FT)/CGAL::exact(Lazy_exact_nt<T>)
#include <CGAL/utils_classes.h>

#ifdef CGAL_LINKED_WITH_TBB
#include <tbb/parallel_reduce.h>
#include <tbb/blocked_range.h>
#endif // CGAL_LINKED_WITH_TBB

#include <boost/container/small_vector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/dynamic_bitset.hpp>

#include <vector>
#include <utility>
#include <algorithm>
#include <unordered_set>

namespace CGAL {

// workaround for area(face_range, tm) overload
template<typename CGAL_NP_TEMPLATE_PARAMETERS_NO_DEFAULT, typename NP>
class GetGeomTraits<CGAL_NP_CLASS, NP>
{
public:
  struct type{};
};

namespace Polygon_mesh_processing {

namespace internal {

inline void rearrange_face_ids(boost::container::small_vector<std::size_t, 4>& ids)
{
  auto min_elem = std::min_element(ids.begin(), ids.end());
  std::rotate(ids.begin(), min_elem, ids.end());
}
}//namespace internal

/**
  * \ingroup PMP_measure_grp
  *
  * computes the length of an edge of a given polygon mesh.
  * The edge is given by one of its halfedges, or the edge itself.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param h one halfedge of the edge whose length is computed
  * @param pmesh the polygon mesh to which `h` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{geom_traits}
  *     \cgalParamDescription{an instance of a geometric traits class}
  *     \cgalParamType{The traits class must provide the nested functor `Compute_squared_distance_3`
  *                    to compute the distance between two points:
  *                    `FT operator()(%Point_3 src1, %Point_3 tgt1)`,
  *                    and a function `Compute_squared_distance_3 compute_squared_distance_3_object()`.}
  *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
  *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
  *   \cgalParamNEnd

  * \cgalNamedParamsEnd
  *
  * @return the length of `h`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `pmesh`.
  *
  * \warning This function involves a square root computation.
  * If `FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @sa `squared_edge_length()`
  * @sa `face_border_length()`
  */
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
edge_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
            const PolygonMesh& pmesh,
            const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Geom_traits;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(is_valid_halfedge_descriptor(h, pmesh));

  typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  return CGAL::approximate_sqrt(gt.compute_squared_distance_3_object()(get(vpm, source(h, pmesh)),
                                                                       get(vpm, target(h, pmesh))));
}

// edge overloads
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
edge_length(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
            const PolygonMesh& pmesh,
            const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(is_valid_edge_descriptor(e, pmesh));

  return edge_length(halfedge(e, pmesh), pmesh, np);
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the squared length of an edge of a given polygon mesh.
  * The edge is given by one of its halfedges, or the edge itself.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param h one halfedge of the edge whose squared length is computed
  * @param pmesh the polygon mesh to which `h` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{geom_traits}
  *     \cgalParamDescription{an instance of a geometric traits class}
  *     \cgalParamType{The traits class must provide the nested functor `Compute_squared_distance_3`
  *                    to compute the distance between two points:
  *                    `FT operator()(%Point_3 src1, %Point_3 tgt1)`,
  *                    and a function `Compute_squared_distance_3 compute_squared_distance_3_object()`.}
  *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
  *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the squared length of `h`. The return type `FT` is a number type either deduced
  * from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `pmesh`.
  *
  * @sa `edge_length()`
  * @sa `face_border_length()`
  */
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
squared_edge_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                    const PolygonMesh& pmesh,
                    const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Geom_traits;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  CGAL_precondition(is_valid_halfedge_descriptor(h, pmesh));

  typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

  return gt.compute_squared_distance_3_object()(get(vpm, source(h, pmesh)),
                                                get(vpm, target(h, pmesh)));
}


// edge overloads
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
squared_edge_length(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
                    const PolygonMesh& pmesh,
                    const NamedParameters& np = parameters::default_values())
{
  CGAL_precondition(is_valid_edge_descriptor(e, pmesh));

  return squared_edge_length(halfedge(e, pmesh), pmesh, np);
}

template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
average_edge_length(const PolygonMesh& pmesh,
                    const NamedParameters& np = parameters::default_values())
{
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;

  const std::size_t n = edges(pmesh).size();
  CGAL_assertion(n > 0);

  typename GT::FT avg_edge_length = 0;
  for (auto e : edges(pmesh))
    avg_edge_length += edge_length(e, pmesh, np);

  avg_edge_length /= static_cast<typename GT::FT>(n);
  return avg_edge_length;
}

#ifdef CGAL_LINKED_WITH_TBB
namespace internal {

template <typename EdgeRange,
          typename PolygonMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
class MinMaxEdgeLength
{
  const EdgeRange& edge_range;
  const PolygonMesh& pmesh;
  const CGAL_NP_CLASS& np;

  using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
  using edge_iterator = typename boost::graph_traits<PolygonMesh>::edge_iterator;

  using Geom_traits = typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;
  using FT = typename Geom_traits::FT;

  static constexpr bool is_random_access =
    std::is_convertible<typename std::iterator_traits<edge_iterator>::iterator_category,
                        std::random_access_iterator_tag>::value;
  std::shared_ptr<std::vector<edge_iterator> > iterators; // to store iterators, if needed

public:
  edge_descriptor min_edge, max_edge;
  FT min_len, max_len;

  MinMaxEdgeLength(const EdgeRange& edge_range,
                    const PolygonMesh& pmesh,
                    const CGAL_NP_CLASS& np)
    : edge_range(edge_range), pmesh(pmesh), np(np),
      min_len(FT(std::numeric_limits<double>::max())),
      max_len(FT(-1))
  {
    if constexpr (!is_random_access)
    {
      // Store iterators for non-random access ranges
      iterators = std::make_shared<std::vector<edge_iterator> >();
      iterators->reserve(std::distance(edge_range.begin(), edge_range.end()));
      for(edge_iterator it = edge_range.begin(); it != edge_range.end(); ++it)
        iterators->push_back(it);
    }
  }

  MinMaxEdgeLength(MinMaxEdgeLength& m, tbb::split)
    : edge_range(m.edge_range), pmesh(m.pmesh), np(m.np),
      iterators(m.iterators),
      min_len(FT(std::numeric_limits<double>::max())),
      max_len(FT(-1))
  { }

  void operator()(const tbb::blocked_range<std::size_t>& r) {
    for(std::size_t i = r.begin(); i != r.end(); ++i) {
      edge_descriptor e;
      if constexpr (is_random_access)
        e = *(edge_range.begin() + i);
      else
        e = *(iterators->at(i));

      FT sq_len = squared_edge_length(e, pmesh, np);
      if(sq_len < min_len) {
        min_len = sq_len;
        min_edge = e;
      }
      if(sq_len > max_len) {
        max_len = sq_len;
        max_edge = e;
      }
    }
  }

  void join(const MinMaxEdgeLength& other) {
    if(other.min_len < min_len) {
      min_len = other.min_len;
      min_edge = other.min_edge;
    }
    if(other.max_len > max_len) {
      max_len = other.max_len;
      max_edge = other.max_edge;
    }
  }
};

} // namespace internal
#endif // CGAL_LINKED_WITH_TBB

/**
  * \ingroup PMP_measure_grp
  *
  * returns the shortest and longest edges of a range of edges of a given polygon mesh,
  * as well as their respective edge lengths.
  *
  * @tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are
  *                        `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
  * @tparam EdgeRange a model of `Range` whose iterator type is `InputIterator` with value type
  *                   `boost::graph_traits<PolygonMesh>::edge_descriptor`.
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param edge_range a range of edges of `pmesh`
  * @param pmesh the polygon mesh in which the longest edge is searched for
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{geom_traits}
  *     \cgalParamDescription{an instance of a geometric traits class}
  *     \cgalParamType{The traits class must provide the nested functor `Compute_squared_distance_3`
  *                    to compute the distance between two points:
  *                    `FT operator()(%Point_3 src1, %Point_3 tgt1)`,
  *                    and a function `Compute_squared_distance_3 compute_squared_distance_3_object()`.}
  *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
  *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the shortest and longest edge in `pmesh`, along with their respective lengths.
  *
  * @warning This function involves a square root computation.
  *
  * @sa `edge_length()`
  * @sa `squared_edge_length()`
  */
template<typename ConcurrencyTag = CGAL::Sequential_tag,
         typename EdgeRange,
         typename PolygonMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::pair<std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT>,
          std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT> >
#else
auto
#endif
minmax_edge_length(const EdgeRange& edge_range,
                   const PolygonMesh& pmesh,
                   const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using edge_descriptor = typename boost::graph_traits<PolygonMesh>::edge_descriptor;
  using edge_iterator = typename boost::graph_traits<PolygonMesh>::edge_iterator;

  using Geom_traits = typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type;
  using FT = typename Geom_traits::FT;

    edge_iterator first = std::cbegin(edge_range), beyond = std::cend(edge_range);
    if(first == beyond)
    {
      return std::make_pair(std::make_pair(edge_descriptor(), FT(0)),
                            std::make_pair(edge_descriptor(), FT(0)));
    }

#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  // parallel
  if constexpr (std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    internal::MinMaxEdgeLength<EdgeRange, PolygonMesh, CGAL_NP_CLASS> reducer(edge_range, pmesh, np);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, edge_range.size()), reducer);

    return std::make_pair(std::make_pair(reducer.min_edge, CGAL::approximate_sqrt(reducer.min_len)),
                          std::make_pair(reducer.max_edge, CGAL::approximate_sqrt(reducer.max_len)));
  }
  else
#endif
  // sequential
  {
    edge_iterator low = first, high = first, eit = first;
    FT sq_lo, sq_hi;
    sq_lo = sq_hi = squared_edge_length(*eit++, pmesh, np);

    for(; eit!=beyond; ++eit)
    {
      const FT sq_l = squared_edge_length(*eit, pmesh, np);

      if(sq_l < sq_lo) {
        low = eit;
        sq_lo = sq_l;
      }
      if(sq_l > sq_hi) {
        high = eit;
        sq_hi = sq_l;
      }
    }

    CGAL_assertion(low != beyond && high != beyond);
    return std::make_pair(std::make_pair(*low, CGAL::approximate_sqrt(sq_lo)),
                          std::make_pair(*high, CGAL::approximate_sqrt(sq_hi)));
  }
}

/*!
 * \ingroup PMP_measure_grp
 * \brief returns the shortest and longest edges of a given polygon mesh,
 * as well as their respective edge lengths.
 * Equivalent to `minmax_edge_length(edges(pmesh), pmesh, np)`
 */
template<typename ConcurrencyTag = CGAL::Sequential_tag,
         typename PolygonMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::pair<std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT>,
          std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT> >
#else
auto
#endif
minmax_edge_length(const PolygonMesh& pmesh,
                   const CGAL_NP_CLASS& np = parameters::default_values())
{
  return minmax_edge_length<ConcurrencyTag>(edges(pmesh), pmesh, np);
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the length of the border polyline that contains a given halfedge.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param h a halfedge of the border polyline whose length is computed
  * @param pmesh the polygon mesh to which `h` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the length of the sequence of border edges of `face(h, pmesh)`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided, or the geometric traits class deduced
  * from the point property map of `pmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @sa `edge_length()`
  */
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
face_border_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                   const PolygonMesh& pmesh,
                   const NamedParameters& np = parameters::default_values())
{
  using FT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT;
  ::CGAL::internal::Evaluate<FT> evaluate;
  FT result = 0;

  for(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor haf : halfedges_around_face(h, pmesh))
  {
    result += edge_length(haf, pmesh, np);
    evaluate(result);
  }

  return result;
}

/**
  * \ingroup PMP_measure_grp
  *
  * finds the longest border of a given triangulated surface and returns
  * a halfedge that is part of this border as well as the length of this border.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param pmesh the polygon mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return a pair composed of two members:
  *   - `first`: a halfedge on the longest border.
  *     The return type `halfedge_descriptor` is a halfedge descriptor. It is
  *     deduced from the graph traits corresponding to the type `PolygonMesh`.
  *     `first` is among the halfedges reported by `extract_boundary_cycles()`.
  *   - `second`: the length of the longest border
  *     The return type `FT` is a number type either deduced from the `geom_traits`
  *     \ref bgl_namedparameters "Named Parameters" if provided,
  *     or the geometric traits class deduced from the point property map of `pmesh`
  *
  * @warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @see `face_border_length()`
  * @see `extract_boundary_cycles()`
  */
template<typename PolygonMesh,
         typename NamedParameters = parameters::Default_named_parameters>
#ifdef DOXYGEN_RUNNING
std::pair<halfedge_descriptor, FT>
#else
std::pair<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
          typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT>
#endif
longest_border(const PolygonMesh& pmesh,
               const NamedParameters& np = parameters::default_values())
{
  typedef typename CGAL::Kernel_traits<
            typename property_map_value<PolygonMesh, CGAL::vertex_point_t>::type>::Kernel::FT  FT;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                       halfedge_descriptor;

  std::vector<halfedge_descriptor> boundary_cycles;
  extract_boundary_cycles(pmesh, std::back_inserter(boundary_cycles));
  halfedge_descriptor result_halfedge = boost::graph_traits<PolygonMesh>::null_halfedge();
  FT result_len = 0;
  for(halfedge_descriptor h : boundary_cycles)
  {
    FT len = face_border_length(h, pmesh, np);

    if(result_len < len)
    {
      result_len = len;
      result_halfedge = h;
    }
  }
  return std::make_pair(result_halfedge, result_len);
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the area of a face of a given triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param f the face whose area is computed
  * @param tmesh the triangulated surface mesh to which `f` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @pre `f != boost::graph_traits<TriangleMesh>::%null_face()`
  *
  * @return the area of `f`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided, or the geometric traits class deduced
  * from the point property map of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @sa `squared_face_area()`
  * @sa `area()`
  */
template<typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
face_area(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
          const TriangleMesh& tmesh,
          const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_face_descriptor(f, tmesh));

  typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  halfedge_descriptor hd = halfedge(f, tmesh);
  halfedge_descriptor nhd = next(hd, tmesh);

  typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  return approximate_sqrt(traits.compute_squared_area_3_object()(get(vpm, source(hd, tmesh)),
                                                                 get(vpm, target(hd, tmesh)),
                                                                 get(vpm, target(nhd, tmesh))));
}


/**
  * \ingroup PMP_measure_grp
  *
  * computes the squared area of a face of a given triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param f the face whose squared area is computed
  * @param tmesh the triangulated surface mesh to which `f` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @pre `f != boost::graph_traits<TriangleMesh>::%null_face()`
  *
  * @return the squared area of `f`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * @sa `face_area()`
  */
template<typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
squared_face_area(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                  const TriangleMesh& tmesh,
                  const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(is_valid_face_descriptor(f, tmesh));

  typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  halfedge_descriptor hd = halfedge(f, tmesh);
  halfedge_descriptor nhd = next(hd, tmesh);

  typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  return traits.compute_squared_area_3_object()(get(vpm, source(hd, tmesh)),
                                                get(vpm, target(hd, tmesh)),
                                                get(vpm, target(nhd, tmesh)));
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the area of a range of faces of a given triangulated surface mesh.
  *
  * @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`.
          Its iterator type is `InputIterator`.
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param face_range the range of faces of whose area is computed
  * @param tmesh the triangulated surface mesh to which the faces of `face_range` belong
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return sum of face areas of `faces`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @sa `face_area()`
  */
template<typename FaceRange,
         typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
area(FaceRange face_range,
     const TriangleMesh& tmesh,
     const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  using FT = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT;
  FT result = 0;
  ::CGAL::internal::Evaluate<FT> evaluate;

  for(face_descriptor f : face_range)
  {
    result += face_area(f, tmesh, np);
    evaluate(result);
  }

  return result;
}

/**
  * \ingroup PMP_measure_grp
  * computes the surface area of a triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `FaceGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tmesh the triangulated surface mesh
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the surface area of `tmesh`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  *
  * @sa `face_area()`
  */
template<typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
area(const TriangleMesh& tmesh,
     const CGAL_NP_CLASS& np = parameters::default_values())
{
  return area(faces(tmesh), tmesh, np);
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the volume of the domain bounded by a closed triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tmesh the closed triangulated surface mesh bounding the volume
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * @pre `tmesh` is closed
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the volume bounded by `tmesh`.
  * The return type `FT` is a number type either deduced from the `geom_traits`
  * \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  */
template<typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
volume(const TriangleMesh& tmesh,
       const CGAL_NP_CLASS& np = parameters::default_values())
{
  CGAL_assertion(is_triangle_mesh(tmesh));
  CGAL_assertion(is_closed(tmesh));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));
  typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::Point_3 origin(0, 0, 0);

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  using FT = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT;
  ::CGAL::internal::Evaluate<FT> evaluate;

  FT volume = 0;
  typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
      CGAL::vertex_point_t>::type>::Kernel::Compute_volume_3 cv3;

  for(face_descriptor f : faces(tmesh))
  {
    volume += cv3(origin,
                  get(vpm, target(halfedge(f, tmesh), tmesh)),
                  get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh)),
                  get(vpm, target(prev(halfedge(f, tmesh), tmesh), tmesh)));
    evaluate(volume);
  }

  return volume;
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the aspect ratio of a face of a given triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param f the face whose aspect ratio is computed
  * @param tmesh the triangulated surface mesh to which `f` belongs
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @pre `f != boost::graph_traits<TriangleMesh>::%null_face()`
  *
  * @return the aspect ratio of `f`. The return type `FT` is a number type
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not support the `sqrt()` operation, the square root computation
  * will be performed approximately.
  */
template<typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT
#endif
face_aspect_ratio(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                  const TriangleMesh& tmesh,
                  const CGAL_NP_CLASS& np = parameters::default_values())
{
  CGAL_precondition(is_valid_face_descriptor(f, tmesh));
  CGAL_precondition(is_triangle(halfedge(f, tmesh), tmesh));

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor           halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type             Geom_traits;
  typedef typename Geom_traits::FT                                                  FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  halfedge_descriptor h = halfedge(f, tmesh);

  Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));

#if 0
  const FT sq_triangle_area = gt.compute_squared_area_3_object()(get(vpm, source(h, tmesh)),
                                                                 get(vpm, target(h, tmesh)),
                                                                 get(vpm, target(next(h, tmesh), tmesh)));
  const FT sq_d12 = gt.compute_squared_distance_2_object()(get(vpm, source(h, tmesh)),
                                                           get(vpm, target(h, tmesh)));
  const FT sq_d13 = gt.compute_squared_distance_2_object()(get(vpm, source(h, tmesh)),
                                                           get(vpm, target(next(h, tmesh), tmesh)));
  const FT sq_d23 = gt.compute_squared_distance_2_object()(get(vpm, target(h, tmesh)),
                                                           get(vpm, target(next(h, tmesh), tmesh)));

  const FT min_sq_d123 = (std::min)(sq_d12, (std::min)(sq_d13, sq_d23));

  const FT aspect_ratio = 4*sq_triangle_area*min_sq_d123 / (sq_d12*sq_d13*sq_d23);
#else // below requires SQRT
  typedef typename Geom_traits::Line_3                                              Line_3;

  FT sq_max_edge_length = gt.compute_squared_distance_3_object()(get(vpm, source(h, tmesh)),
                                                                 get(vpm, target(h, tmesh)));
  FT sq_min_alt = gt.compute_squared_distance_3_object()(get(vpm, target(next(h, tmesh), tmesh)),
                                                         Line_3(get(vpm, source(h, tmesh)),
                                                                get(vpm, target(h, tmesh))));
  h = next(h, tmesh);

  for(int i=1; i<3; ++i)
  {
    FT sq_edge_length = gt.compute_squared_distance_3_object()(get(vpm, source(h, tmesh)),
                                                               get(vpm, target(h, tmesh)));
    FT sq_alt = gt.compute_squared_distance_3_object()(get(vpm, target(next(h, tmesh), tmesh)),
                                                       Line_3(get(vpm, source(h, tmesh)),
                                                              get(vpm, target(h, tmesh))));

    if(sq_alt < sq_min_alt)
      sq_min_alt = sq_alt;
    if(sq_edge_length > sq_max_edge_length)
      sq_max_edge_length = sq_edge_length;

    h = next(h, tmesh);
  }

  CGAL_assertion(sq_min_alt > 0);
  const FT aspect_ratio = CGAL::approximate_sqrt(sq_max_edge_length / sq_min_alt);
#endif

  return aspect_ratio;
}

/**
  * \ingroup PMP_measure_grp
  *
  * computes the centroid of a volume bounded by a closed triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `FaceListGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param tmesh the closed triangulated surface mesh bounding the volume
  * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * @pre `tmesh` is closed
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`}
  *   \cgalParamNEnd
  *
 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
  * @return the centroid of the domain bounded by `tmesh`.
  */
template<typename TriangleMesh, typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
Point_3
#else
typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::Point_3
#endif
centroid(const TriangleMesh& tmesh,
         const CGAL_NP_CLASS& np = parameters::default_values())
{
  // See: http://www2.imperial.ac.uk/~rn/centroid.pdf

  CGAL_assertion(is_triangle_mesh(tmesh));
  CGAL_assertion(is_closed(tmesh));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type     Kernel;
  Kernel k = choose_parameter<Kernel>(get_parameter(np, internal_np::geom_traits));

  typedef typename Kernel::FT                                           FT;
  typedef typename boost::property_traits<Vpm>::reference               Point_3_ref;
  typedef typename Kernel::Vector_3                                     Vector_3;

  typedef typename Kernel::Construct_translated_point_3                 Construct_translated_point_3;
  typedef typename Kernel::Construct_vector_3                           Construct_vector_3;
  typedef typename Kernel::Construct_normal_3                           Construct_normal_3;
  typedef typename Kernel::Compute_scalar_product_3                     Scalar_product;
  typedef typename Kernel::Construct_scaled_vector_3                    Scale;
  typedef typename Kernel::Construct_sum_of_vectors_3                   Sum;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor   face_descriptor;

  FT volume = 0;

  Vector_3 centroid(NULL_VECTOR);

  Construct_translated_point_3 point = k.construct_translated_point_3_object();
  Construct_vector_3 vector = k.construct_vector_3_object();
  Construct_normal_3 normal = k.construct_normal_3_object();
  Scalar_product scalar_product = k.compute_scalar_product_3_object();
  Scale scale = k.construct_scaled_vector_3_object();
  Sum sum = k.construct_sum_of_vectors_3_object();

  for(face_descriptor fd : faces(tmesh))
  {
    const Point_3_ref p = get(vpm, target(halfedge(fd, tmesh), tmesh));
    const Point_3_ref q = get(vpm, target(next(halfedge(fd, tmesh), tmesh), tmesh));
    const Point_3_ref r = get(vpm, target(prev(halfedge(fd, tmesh), tmesh), tmesh));
    Vector_3 vp = vector(ORIGIN, p),
             vq = vector(ORIGIN, q),
             vr = vector(ORIGIN, r);
    Vector_3 n = normal(p, q, r);
    volume += (scalar_product(n,vp))/FT(6);
    n = scale(n, FT(1)/FT(24));

    Vector_3 v2 = sum(vp, vq);
    Vector_3 v3 = Vector_3(square(v2.x()), square(v2.y()), square(v2.z()));
    v2 = sum(vq, vr);
    v3 = sum(v3, Vector_3(square(v2.x()), square(v2.y()), square(v2.z())));
    v2 = sum(vp, vr);
    v3 = sum(v3, Vector_3(square(v2.x()), square(v2.y()), square(v2.z())));

    centroid = sum(centroid, Vector_3(n.x() * v3.x(), n.y() * v3.y(), n.z() * v3.z()));
  }

  centroid = scale(centroid, FT(1)/(FT(2)*volume));
  return point(ORIGIN, centroid);
}

/**
  * \ingroup PMP_measure_grp
  *
  * identifies faces only present in `m1` and `m2` as well as the faces present
  * in both polygon meshes. Two faces are matching if they have the same
  * orientation and the same points.
  *
  * @tparam PolygonMesh1 a model of `HalfedgeListGraph` and `FaceListGraph`
  * @tparam PolygonMesh2 a model of `HalfedgeListGraph` and `FaceListGraph`
  * @tparam FaceOutputIterator1 model of `OutputIterator`
  *   holding `boost::graph_traits<PolygonMesh1>::%face_descriptor`.
  * @tparam FaceOutputIterator2 model of `OutputIterator`
  *   holding `boost::graph_traits<PolygonMesh2>::%face_descriptor`.
  * @tparam FacePairOutputIterator model of `OutputIterator`
  *   holding `std::pair<boost::graph_traits<PolygonMesh1>::%face_descriptor,
  *   boost::graph_traits<PolygonMesh2>::%face_descriptor`.
  *
  * @tparam NamedParameters1 a sequence of \ref bgl_namedparameters "Named Parameters"
  * @tparam NamedParameters2 a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param m1 the first polygon mesh
  * @param m2 the second polygon mesh
  * @param common output iterator collecting the faces that are common to both meshes
  * @param m1_only output iterator collecting the faces that are only in `m1`
  * @param m2_only output iterator collecting the faces that are only in `m2`
  * @param np1 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  * @param np2 an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
  *
  * \cgalNamedParamsBegin
  *   \cgalParamNBegin{vertex_point_map}
  *     \cgalParamDescription{a property map associating points to the vertices of `m1`}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh1>::%vertex_descriptor`
  *                    as key type and `%Point_3` as value type. `%Point_3` must be `LessThanComparable`.}
  *     \cgalParamDefault{`boost::get(CGAL::vertex_point, m1)`}
  *     \cgalParamExtra{The same holds for `m2` and `PolygonMesh2` and the point type must be the same for both meshes.}
  *   \cgalParamNEnd
  *
  *   \cgalParamNBegin{vertex_index_map}
  *     \cgalParamDescription{a property map associating to each vertex of `m1` a unique index between `0` and `num_vertices(m1) - 1`, and similarly for `m2`.}
  *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<Graph>::%vertex_descriptor`
  *                    as key type and `std::size_t` as value type}
  *     \cgalParamDefault{an automatically indexed internal map}
  *     \cgalParamExtra{If this parameter is not passed, internal machinery will create and initialize
  *                     a face index property map, either using the internal property map if it exists
  *                     or using an external map. The latter might result in  - slightly - worsened performance
  *                     in case of non-constant complexity for index access. The same holds for `m2` and `PolygonMesh2`.}
  *   \cgalParamNEnd
  * \cgalNamedParamsEnd
  *
 */
template< typename PolygonMesh1,
          typename PolygonMesh2,
          typename FacePairOutputIterator,
          typename FaceOutputIterator1,
          typename FaceOutputIterator2,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters >
void match_faces(const PolygonMesh1& m1,
                 const PolygonMesh2& m2,
                 FacePairOutputIterator common,
                 FaceOutputIterator1 m1_only,
                 FaceOutputIterator2 m2_only,
                 const NamedParameters1& np1 = parameters::default_values(),
                 const NamedParameters2& np2 = parameters::default_values())
{
  typedef typename GetVertexPointMap<PolygonMesh1, NamedParameters1>::const_type            VPMap1;
  typedef typename GetVertexPointMap<PolygonMesh2, NamedParameters2>::const_type            VPMap2;
  typedef typename GetInitializedVertexIndexMap<PolygonMesh1, NamedParameters1>::const_type VIMap1;
  typedef typename GetInitializedVertexIndexMap<PolygonMesh2, NamedParameters2>::const_type VIMap2;
  typedef typename boost::property_traits<VPMap2>::value_type                               Point_3;
  typedef typename boost::graph_traits<PolygonMesh1>::face_descriptor                       face_descriptor_1;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const VPMap1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                       get_const_property_map(vertex_point, m1));
  const VPMap2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                       get_const_property_map(vertex_point, m2));
  static_assert(std::is_same<typename boost::property_traits<VPMap1>::value_type,
                             typename boost::property_traits<VPMap2>::value_type>::value,
                            "Both vertex point maps must have the same point type.");

  const VIMap1 vim1 = get_initialized_vertex_index_map(m1, np1);
  const VIMap2 vim2 = get_initialized_vertex_index_map(m2, np2);

  std::map<Point_3, std::size_t> point_id_map;

  std::vector<std::size_t> m1_vertex_id(num_vertices(m1), -1);
  std::vector<std::size_t> m2_vertex_id(num_vertices(m2), -1);
  boost::dynamic_bitset<> shared_vertices(m1_vertex_id.size() + m2_vertex_id.size());

  //iterate both meshes to set ids of all points, and set vertex/point_id maps.
  std::size_t id = 0;
  for(auto v : vertices(m1))
  {
    const typename boost::property_traits<VPMap1>::reference p = get(vpm1, v);
    auto res = point_id_map.emplace(p, id);
    if(res.second)
      ++id;
    m1_vertex_id[get(vim1, v)] = res.first->second;
  }
  for(auto v : vertices(m2))
  {
    const typename boost::property_traits<VPMap2>::reference p = get(vpm2, v);
    auto res = point_id_map.emplace(p, id);
    if(res.second)
      ++id;
    else
      shared_vertices.set(res.first->second);
    m2_vertex_id[get(vim2, v)] = res.first->second;
  }

  //fill a set with the "faces point-ids" of m1 and then iterate faces of m2 to compare.
  std::map<boost::container::small_vector<std::size_t, 4>, face_descriptor_1> m1_faces_map;
  for(auto f : faces(m1))
  {
    bool all_shared = true;
    boost::container::small_vector<std::size_t, 4> ids;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m1), m1))
    {
      std::size_t vid = m1_vertex_id[get(vim1, v)];
      ids.push_back(vid);
      if(!shared_vertices.test(vid))
      {
        all_shared = false;
        break;
      }
    }
    if(all_shared)
    {
      internal::rearrange_face_ids(ids);
      m1_faces_map.emplace(ids, f);
    }
    else
      *m1_only++ = f;
  }
  for(auto f : faces(m2))
  {
    boost::container::small_vector<std::size_t, 4> ids;
    bool all_shared = true;
    for(auto v : CGAL::vertices_around_face(halfedge(f, m2), m2))
    {
      std::size_t vid = m2_vertex_id[get(vim2, v)];
      ids.push_back(vid);
      if(!shared_vertices.test(vid))
      {
        all_shared = false;
        break;
      }
    }
    if(all_shared)
    {
      internal::rearrange_face_ids(ids);
      auto it = m1_faces_map.find(ids);
      if(it != m1_faces_map.end())
      {
        *common++ = std::make_pair(it->second, f);
        m1_faces_map.erase(it);
      }
      else
      {
        *m2_only++ = f;
      }
    }
    else
      *m2_only++ = f;
  }
  //all shared faces have been removed from the map, so all that remains must go in m1_only
  for(const auto& it : m1_faces_map)
  {
    *m1_only++ = it.second;
  }
}

#ifdef CGAL_LINKED_WITH_TBB
namespace internal {

template <typename EdgeRange,
          typename TriangleMesh,
          typename CGAL_NP_TEMPLATE_PARAMETERS>
class MinMaxDihedralAngle
{
  const EdgeRange& edge_range;
  const TriangleMesh& tmesh;
  const CGAL_NP_CLASS& np;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using edge_iterator = typename boost::graph_traits<TriangleMesh>::edge_iterator;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type;
  using FT = typename Geom_traits::FT;

  static constexpr bool is_random_access =
    std::is_convertible<typename std::iterator_traits<edge_iterator>::iterator_category,
                        std::random_access_iterator_tag>::value;
  std::shared_ptr<std::vector<edge_iterator> > iterators; // to store iterators, if needed

public:
  edge_descriptor min_edge, max_edge;
  FT min_angle, max_angle;

  MinMaxDihedralAngle(const EdgeRange& edge_range,
                       const TriangleMesh& tmesh,
                       const CGAL_NP_CLASS& np)
    : edge_range(edge_range), tmesh(tmesh), np(np),
      min_angle(200), max_angle(-200)
  {
    if constexpr (!is_random_access)
    {
      // Store iterators for non-random access ranges
      iterators = std::make_shared<std::vector<edge_iterator> >();
      iterators->reserve(std::distance(edge_range.begin(), edge_range.end()));
      for(edge_iterator it = edge_range.begin(); it != edge_range.end(); ++it)
        iterators->push_back(it);
    }
  }

  MinMaxDihedralAngle(MinMaxDihedralAngle& m, tbb::split)
    : edge_range(m.edge_range), tmesh(m.tmesh), np(m.np),
      iterators(m.iterators),
      min_angle(200), max_angle(-200)
  { }

  void operator()(const tbb::blocked_range<std::size_t>& r)
  {
    using parameters::choose_parameter;
    using parameters::get_parameter;

    using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;

    typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, tmesh));

    Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
    typename Geom_traits::Compute_approximate_dihedral_angle_3 approx_dh =
      gt.compute_approximate_dihedral_angle_3_object();

    for(std::size_t i = r.begin(); i != r.end(); ++i)
    {
      edge_descriptor e;
      if constexpr (is_random_access)
        e = *(edge_range.begin() + i);
      else
        e = *(iterators->at(i));

      if(is_border(e, tmesh))
        continue;

      const halfedge_descriptor h = halfedge(e, tmesh);
      CGAL_assertion(!is_degenerate_triangle_face(h, tmesh));

      const vertex_descriptor p = source(h,tmesh);
      const vertex_descriptor q = target(h,tmesh);
      const vertex_descriptor r = target(next(h,tmesh),tmesh);
      const vertex_descriptor s = target(next(opposite(h,tmesh),tmesh),tmesh);

      const FT angle = approx_dh(get(vpm,p), get(vpm,q), get(vpm,r), get(vpm,s));
      if(angle < min_angle) {
        min_angle = angle;
        min_edge = e;
      }
      if(angle > max_angle) {
        max_angle = angle;
        max_edge = e;
      }
    }
  }

  void join(const MinMaxDihedralAngle& other)
  {
    if(other.min_angle < min_angle) {
      min_angle = other.min_angle;
      min_edge = other.min_edge;
    }
    if(other.max_angle > max_angle) {
      max_angle = other.max_angle;
      max_edge = other.max_edge;
    }
  }
};

} // namespace internal
#endif // CGAL_LINKED_WITH_TBB

/*!
 * \ingroup PMP_measure_grp
 *
 * \brief computes the minimum and maximum dihedral angles of a range of edges of a given polygon mesh.
 *
 * \tparam ConcurrencyTag enables sequential versus parallel algorithm. Possible values are
 *                        `Sequential_tag`, `Parallel_tag`, and `Parallel_if_available_tag`.
 * \tparam EdgeRange a model of `Range` whhose iterator type is `InputIterator` with value type
 *                   `boost::graph_traits<PolygonMesh>::edge_descriptor`.
 * \tparam TriangleMesh a model of `HalfedgeListGraph`
 * \tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
 *
 * \param edge_range a range of edges of `tmesh`
 * \param tmesh the polygon mesh
 * \param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below.
 *
 * \cgalNamedParamsBegin
 *
 *   \cgalParamNBegin{vertex_point_map}
 *     \cgalParamDescription{a property map associating points to the vertices of `tmesh`.}
 *     \cgalParamType{a class model of `ReadablePropertyMap` with
 *                    `boost::graph_traits<TriangleMesh>::%vertex_descriptor`
 *                    as key type and `geom_traits::Point_3` as value type.}
 *     \cgalParamDefault{`boost::get(CGAL::vertex_point, tmesh)`.}
 *     \cgalParamExtra{If this parameter is omitted, an internal property map for
 *                     `CGAL::vertex_point_t` must be available in `TriangleMesh`.}
 *   \cgalParamNEnd

 *   \cgalParamNBegin{geom_traits}
 *     \cgalParamDescription{an instance of a geometric traits class}
 *     \cgalParamType{a class model of `Kernel`}
 *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
 *     \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
 *   \cgalParamNEnd
 * \cgalNamedParamsEnd
 *
 * \pre `CGAL::is_triangle_mesh(tmesh)`
 * \pre There are no degenerate faces in `tmesh`.
 *
 * \see `detect_sharp_edges()`
 * \see `sharp_edges_segmentation()`
 */
template<typename ConcurrencyTag = CGAL::Sequential_tag,
         typename EdgeRange,
         typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::pair<std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT>,
          std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT> >
#else
auto
#endif
minmax_dihedral_angle(const EdgeRange& edge_range,
                      const TriangleMesh& tmesh,
                      const CGAL_NP_CLASS& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
  using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
  using edge_iterator = typename boost::graph_traits<TriangleMesh>::edge_iterator;

  using Geom_traits = typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type;
  using FT = typename Geom_traits::FT;

  edge_iterator first = std::cbegin(edge_range), beyond = std::cend(edge_range);
  if(first == beyond)
  {
    return std::make_pair(std::make_pair(edge_descriptor(), FT(200)),
                          std::make_pair(edge_descriptor(), FT(-200)));
  }

#if !defined(CGAL_LINKED_WITH_TBB)
  static_assert (!(std::is_convertible<ConcurrencyTag, Parallel_tag>::value),
                             "Parallel_tag is enabled but TBB is unavailable.");
#else
  // parallel
  if constexpr (std::is_convertible<ConcurrencyTag, Parallel_tag>::value)
  {
    internal::MinMaxDihedralAngle<EdgeRange, TriangleMesh, CGAL_NP_CLASS> reducer(edge_range, tmesh, np);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, edge_range.size()), reducer);

    return std::make_pair(std::make_pair(reducer.min_edge, reducer.min_angle),
                          std::make_pair(reducer.max_edge, reducer.max_angle));
  }
  else
#endif
  // sequential
  {
    edge_descriptor min_edge, max_edge;
    FT min_angle(200), max_angle(-200);

    typename GetVertexPointMap<TriangleMesh, CGAL_NP_CLASS>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, tmesh));

    Geom_traits gt = choose_parameter<Geom_traits>(get_parameter(np, internal_np::geom_traits));
    typename Geom_traits::Compute_approximate_dihedral_angle_3 approx_dh =
      gt.compute_approximate_dihedral_angle_3_object();

    for(edge_iterator eit = first; eit != beyond; ++eit)
    {
      if(is_border(*eit, tmesh))
        continue;

      const halfedge_descriptor h = halfedge(*eit, tmesh);
      CGAL_assertion(!is_degenerate_triangle_face(h, tmesh));

       const vertex_descriptor p = source(h,tmesh);
       const vertex_descriptor q = target(h,tmesh);
       const vertex_descriptor r = target(next(h,tmesh),tmesh);
       const vertex_descriptor s = target(next(opposite(h,tmesh),tmesh),tmesh);

       const FT angle = approx_dh(get(vpm,p), get(vpm,q), get(vpm,r), get(vpm,s));
      if(angle < min_angle) {
        min_angle = angle;
        min_edge = *eit;
      }
      if(angle > max_angle) {
        max_angle = angle;
        max_edge = *eit;
      }
    }

    return std::make_pair(std::make_pair(min_edge, min_angle),
                          std::make_pair(max_edge, max_angle));
  }
}

/*!
 * \ingroup PMP_measure_grp
 * computes the minimum and maximum dihedral angles of a given triangle mesh.
 * Equivalent to `minmax_dihedral_angle(edges(tmesh), tmesh, np)`
 */
template<typename ConcurrencyTag = CGAL::Sequential_tag,
         typename TriangleMesh,
         typename CGAL_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
std::pair<std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT>,
          std::pair<boost::graph_traits<PolygonMesh>::edge_descriptor`, FT> >
#else
auto
#endif
minmax_dihedral_angle(const TriangleMesh& tmesh,
                      const CGAL_NP_CLASS& np = parameters::default_values())
{
  return minmax_dihedral_angle<ConcurrencyTag>(edges(tmesh), tmesh, np);
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif // CGAL_POLYGON_MESH_PROCESSING_MEASURE_H
