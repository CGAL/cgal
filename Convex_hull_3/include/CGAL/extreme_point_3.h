// Copyright (c) 2022 INRIA (France) and GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Léo Valque
//

#ifndef CGAL_EXTREME_POINT_3_H
#define CGAL_EXTREME_POINT_3_H

#include <CGAL/license/Convex_hull_3.h>

#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Kernel_23/internal/Has_boolean_tags.h>

#include <CGAL/Container_helper.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/property_map.h>

#include <boost/range/value_type.hpp>
#include <boost/range/reference.hpp>

#include <boost/graph/adjacency_list.hpp>

#include <vector>

namespace CGAL {

/**
* \ingroup PkgConvexHull3Functions
*
* computes the furthest point of the range along the direction.
*
* @tparam PointRange is a model of `ConstRange`. The value type of its iterator is the key type of the named parameter `point_map`.
* @tparam Direction_3 is a model of CGAL::Direction_3.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param r the range of points
* @param dir the direction
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{point_map}
*     \cgalParamDescription{a property map associating points to the elements of `r`}
*     \cgalParamType{a model of `ReadablePropertyMap` whose value types are the same for `np1` and `np2`}
*     \cgalParamDefault{`CGAL::Identity_property_map`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from `Direction_3`, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits_converter}
*     \cgalParamDescription{A Converter from the point type of `point_map` to the point type of `geom_traits`}
*     \cgalParamType{a class model of `NT_Converter`}
*     \cgalParamDefault{a \cgal `Cartesian_converter` deduced from ` point_map` and `geom_traits`, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \return an instance of the key type of `point_map` parameter.
*/
template <class Range, class Direction_3, class NamedParameters>
#if DOXYGEN_RUNNING
Point_3
#else
typename Point_set_processing_3_np_helper<Range, NamedParameters>::Const_point_map::key_type
#endif
// typename Kernel_traits<Vector_3>::Kernel::Point_3
extreme_point_3(const Range& r, const Direction_3 &dir, const NamedParameters &np) {
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;
  using NP_helper= Point_set_processing_3_np_helper<Range, NamedParameters>;

  using PointMap= typename NP_helper::Const_point_map;
  PointMap point_map = NP_helper::get_const_point_map(r, np);

  using Point_GT = typename NP_helper::Geom_traits;
  using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
  using GT=typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParameters,
      Default_GT
    > ::type;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using Default_geom_traits_converter = Cartesian_converter<Point_GT, GT>;
  using GTC=typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_converter_t,
      NamedParameters,
      Default_geom_traits_converter
    > ::type;
  GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));

  auto vector_3 = gt.construct_vector_3_object();
  auto csp = gt.compare_scalar_product_3_object();

  typename Range::const_iterator argmax=r.begin();
  Vector_3 vec_max = vector_3(ORIGIN, converter(get(point_map, *argmax)));
  for(typename Range::const_iterator it=++r.begin(); it!=r.end(); ++it){
    Vector_3 vec = vector_3(ORIGIN, converter(get(point_map, *it)));
    if(csp(vec_max, dir.vector(), vec, dir.vector())==SMALLER){
      vec_max=vec;
      argmax=it;
    }
  }
  return *argmax;
}

/**
* \ingroup PkgConvexHull3Functions
*
* computes the furthest point of the convex graph along the direction `dir` and returns the associated vertex.
*
* @tparam Graph is a model of `VertexListGraph` and `AdjacencyGraph`.
* @tparam Direction_3 is a model of `Kernel::Direction_3`.
* @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
*
* @param g the convex graph
* @param dir the direction
* @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
*
* \pre The input graph must represent a convex object to guarantee a correct answer.
*
* \cgalNamedParamsBegin
*   \cgalParamNBegin{vertex_point_map}
*     \cgalParamDescription{a property map associating points to the vertices of `g`}
*     \cgalParamType{a model of `ReadablePropertyMap`}
*     \cgalParamDefault{boost::get(CGAL::vertex_point, g)}
*     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t` must be available in `Graph`.}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits}
*     \cgalParamDescription{An instance of a geometric traits class}
*     \cgalParamType{a class model of `Kernel`}
*     \cgalParamDefault{a \cgal kernel deduced from `Direction_3`, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
*   \cgalParamNBegin{geom_traits_converter}
*     \cgalParamDescription{A converter from the point type of `vertex_point_map` to the point type of `geom_traits`}
*     \cgalParamType{a class model of `NT_Converter`}
*     \cgalParamDefault{a \cgal `Cartesian_converter` deduced from ` vertex_point_map`  and `geom_traits`, using `CGAL::Kernel_traits`}
*   \cgalParamNEnd
* \cgalNamedParamsEnd
*
* \return a `boost::graph_traits<Graph>::vertex_descriptor`
*/
template <class Graph, class NamedParameters, class Direction_3>
#if DOXYGEN_RUNNING
vertex_descriptor
#else
typename boost::graph_traits<Graph>::vertex_descriptor
#endif
extreme_vertex_3(const Graph& g, const Direction_3 &dir, const NamedParameters &np) {
  using CGAL::parameters::choose_parameter;
  using CGAL::parameters::get_parameter;

  using vertex_descriptor= typename boost::graph_traits<Graph>::vertex_descriptor;

  using GetVertexPointMap = GetVertexPointMap<Graph, NamedParameters>;
  using VPM = typename GetVertexPointMap::const_type;
  VPM point_map = GetVertexPointMap::get_const_map(np, g);

  using GetGeomTraits = GetGeomTraits<Graph, NamedParameters>;
  using GraphGT= typename GetGeomTraits::type;

  using Default_GT = typename Kernel_traits<Direction_3>::Kernel;
  using GT=typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_t,
      NamedParameters,
      Default_GT
    > ::type;
  GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  using Default_geom_traits_converter = Cartesian_converter<GraphGT, GT>;
  using GTC=typename internal_np::Lookup_named_param_def <
      internal_np::geom_traits_converter_t,
      NamedParameters,
      Default_geom_traits_converter
    > ::type;
  GTC converter = choose_parameter<GTC>(get_parameter(np, internal_np::geom_traits_converter));

  auto vector_3 = gt.construct_vector_3_object();
  auto csp = gt.compare_scalar_product_3_object();

  // If the number of vertices is small, simply test all vertices
  if(vertices(g).size()<20)
    return extreme_point_3(vertices(g), dir, parameters::point_map(point_map));

  //Walks on the mesh to find a local maximun
  vertex_descriptor argmax = *vertices(g).begin();
  Vector_3 vec_max = vector_3(ORIGIN, converter(get(point_map,argmax)));
  bool is_local_max;
  do{
    is_local_max=true;
    for(auto v: vertices_around_target(argmax, g)){
      Vector_3 vec = vector_3(ORIGIN, converter(get(point_map, v)));
      if(csp(vec_max, dir.vector(), vec, dir.vector())==SMALLER){
        vec_max = vec;
        argmax = v;
        is_local_max=false; // repeat with the new vertex
        break;
      }
    }
  }while(!is_local_max);
  // Since convex, local maximum is a global maximum
  return argmax;
}

// #endif

} // CGAL namespace

#endif // CGAL_EXTREME_POINT_3_H
