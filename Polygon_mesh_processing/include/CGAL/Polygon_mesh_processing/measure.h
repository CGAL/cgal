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

#include <CGAL/disable_warnings.h>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <CGAL/Polygon_mesh_processing/border.h>

#include <CGAL/Lazy.h> // needed for CGAL::exact(FT)/CGAL::exact(Lazy_exact_nt<T>)

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
 *     \cgalParamType{a class model of `Kernel`}
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
 *     \cgalParamType{a class model of `Kernel`}
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
  typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT result = 0;

  for(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor haf : halfedges_around_face(h, pmesh))
  {
    result += edge_length(haf, pmesh, np);
    exact(result);
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

  typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT result = 0;
  for(face_descriptor f : face_range)
  {
    result += face_area(f, tmesh, np);
    exact(result);
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

  typename GetGeomTraits<TriangleMesh, CGAL_NP_CLASS>::type::FT volume = 0;
  typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
      CGAL::vertex_point_t>::type>::Kernel::Compute_volume_3 cv3;

  for(face_descriptor f : faces(tmesh))
  {
    volume += cv3(origin,
                  get(vpm, target(halfedge(f, tmesh), tmesh)),
                  get(vpm, target(next(halfedge(f, tmesh), tmesh), tmesh)),
                  get(vpm, target(prev(halfedge(f, tmesh), tmesh), tmesh)));
    exact(volume);
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

} // namespace Polygon_mesh_processing
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_MEASURE_H
