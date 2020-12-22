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
#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/Kernel/global_functions_3.h>

#include <CGAL/Lazy.h> // needed for CGAL::exact(FT)/CGAL::exact(Lazy_exact_nt<T>)

#include <boost/unordered_set.hpp>
#include <boost/graph/graph_traits.hpp>

#include <utility>

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {

// workaround for area(face_range, tm) overload
template<typename CGAL_PMP_NP_TEMPLATE_PARAMETERS, typename NP>
class GetGeomTraits<CGAL_PMP_NP_CLASS, NP>
{
public:
  struct type{};
};

namespace Polygon_mesh_processing {

/**
  * \ingroup measure_grp
  * computes the length of an edge of a given polygon mesh.
  * The edge is given by one of its halfedges, or the edge itself.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param h one halfedge of the edge to compute the length
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
  * @return the length of `h`. The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `pmesh`.
  *
  * \warning This function involves a square root computation.
  * If `FT` does not have a `sqrt()` operation, the square root computation
  * will be done approximately.
  *
  * @sa `face_border_length()`
  */
template<typename PolygonMesh,
         typename NamedParameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
edge_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
            const PolygonMesh& pmesh,
            const NamedParameters& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, pmesh));

  return CGAL::approximate_sqrt(CGAL::squared_distance(get(vpm, source(h, pmesh)),
                                                       get(vpm, target(h, pmesh))));
}

template<typename PolygonMesh>
typename CGAL::Kernel_traits<typename property_map_value<PolygonMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
edge_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
            const PolygonMesh& pmesh)
{
  return edge_length(h, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}
// edge overloads
template<typename PolygonMesh,
         typename NamedParameters>
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
edge_length(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
            const PolygonMesh& pmesh,
            const NamedParameters& np)
{
  return edge_length(halfedge(e, pmesh), pmesh, np);
}

template<typename PolygonMesh>
typename CGAL::Kernel_traits<typename property_map_value<PolygonMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
edge_length(typename boost::graph_traits<PolygonMesh>::edge_descriptor e,
            const PolygonMesh& pmesh)
{
  return edge_length(halfedge(e, pmesh), pmesh);
}

/**
  * \ingroup measure_grp
  * computes the length of the border polyline
  * that contains a given halfedge.
  *
  * @tparam PolygonMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param h a halfedge of the border polyline of which the length is computed
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
  * The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `pmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
  * will be done approximately.
  *
  * @sa `edge_length()`
  */
template<typename PolygonMesh,
         typename NamedParameters>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT
#endif
face_border_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                   const PolygonMesh& pmesh,
                   const NamedParameters& np)
{
  typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT result = 0;

  for(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor haf : halfedges_around_face(h, pmesh))
  {
    result += edge_length(haf, pmesh, np);
    exact(result);
  }

  return result;
}

template<typename PolygonMesh>
typename CGAL::Kernel_traits<typename property_map_value<PolygonMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
face_border_length(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor h,
                   const PolygonMesh& pmesh)
{
  return face_border_length(h, pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * finds the longest border of a given triangulated surface and returns
  * a halfedge that is part of this border and the length of this border.
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
  *   - `second`: the length of the longest border
  *     The return type `FT` is a number type. It is
  *     either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  *     or the geometric traits class deduced from the point property map
  *     of `pmesh`
  *
  */
template<typename PolygonMesh,
         typename NamedParameters>
#ifdef DOXYGEN_RUNNING
std::pair<halfedge_descriptor, FT>
#else
std::pair<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
typename GetGeomTraits<PolygonMesh, NamedParameters>::type::FT>
#endif
longest_border(const PolygonMesh& pmesh,
               const NamedParameters& np)
{
  typedef typename CGAL::Kernel_traits<
            typename property_map_value<PolygonMesh, CGAL::vertex_point_t>::type>::Kernel::FT  FT;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor                       halfedge_descriptor;

  boost::unordered_set<halfedge_descriptor> visited;
  halfedge_descriptor result_halfedge = boost::graph_traits<PolygonMesh>::null_halfedge();
  FT result_len = 0;
  for(halfedge_descriptor h : halfedges(pmesh))
  {
    if(visited.find(h)== visited.end())
    {
      if(is_border(h, pmesh))
      {
        FT len = 0;
        for(halfedge_descriptor haf : halfedges_around_face(h, pmesh))
        {
          len += edge_length(haf, pmesh, np);
          visited.insert(haf);
        }

        if(result_len < len)
        {
          result_len = len;
          result_halfedge = h;
        }
      }
    }
  }
  return std::make_pair(result_halfedge, result_len);
}

template<typename PolygonMesh>
std::pair<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor,
typename CGAL::Kernel_traits<typename property_map_value<PolygonMesh,
CGAL::vertex_point_t>::type>::Kernel::FT>
longest_border(const PolygonMesh& pmesh)
{
  return longest_border(pmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the area of a face of a given
  * triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param f the face of which the area is computed
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
  * The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `tmesh`.
  *
  * @sa `area()`
  */
template<typename TriangleMesh,
         typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT
#endif
face_area(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
          const TriangleMesh& tmesh,
          const CGAL_PMP_NP_CLASS& np)
{
  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  CGAL_precondition(boost::graph_traits<TriangleMesh>::null_face() != f);

  typename GetVertexPointMap<TriangleMesh, CGAL_PMP_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  halfedge_descriptor hd = halfedge(f, tmesh);
  halfedge_descriptor nhd = next(hd, tmesh);

  typedef typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type GT;
  GT traits = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));

  return approximate_sqrt(traits.compute_squared_area_3_object()(get(vpm, source(hd, tmesh)),
                                                                 get(vpm, target(hd, tmesh)),
                                                                 get(vpm, target(nhd, tmesh))));
}

template<typename TriangleMesh>
typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
face_area(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
          const TriangleMesh& tmesh)
{
  return face_area(f, tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the area of a range of faces of a given
  * triangulated surface mesh.
  *
  * @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`.
          Its iterator type is `InputIterator`.
  * @tparam TriangleMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param face_range the range of faces of which the area is computed
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
  * The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
  * will be done approximately.
  *
  * @sa `face_area()`
  */
template<typename FaceRange,
         typename TriangleMesh,
         typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT
#endif
area(FaceRange face_range,
     const TriangleMesh& tmesh,
     const CGAL_PMP_NP_CLASS& np)
{
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT result = 0;
  for(face_descriptor f : face_range)
  {
    result += face_area(f, tmesh, np);
    exact(result);
  }

  return result;
}

template<typename FaceRange, typename TriangleMesh>
typename GetGeomTraits<TriangleMesh>::type::FT
area(FaceRange face_range, const TriangleMesh& tmesh)
{
  return area(face_range, tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the surface area of a triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `HalfedgeGraph`
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
  * The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `tmesh`.
  *
  * \warning This function involves a square root computation.
  * If `Kernel::FT` does not have a `sqrt()` operation, the square root computation
  * will be done approximately.
  *
  * @sa `face_area()`
  */
template<typename TriangleMesh,
         typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT
#endif
area(const TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  return area(faces(tmesh), tmesh, np);
}

template<typename TriangleMesh>
typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
area(const TriangleMesh& tmesh)
{
  return area(faces(tmesh), tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the volume of the domain bounded by
  * a closed triangulated surface mesh.
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
  * The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map
  * of `tmesh`.
  */
template<typename TriangleMesh,
         typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT
#endif
volume(const TriangleMesh& tmesh,
       const CGAL_PMP_NP_CLASS& np)
{
  CGAL_assertion(is_triangle_mesh(tmesh));
  CGAL_assertion(is_closed(tmesh));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename GetVertexPointMap<TriangleMesh, CGAL_PMP_NP_CLASS>::const_type
      vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));
  typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::Point_3 origin(0, 0, 0);

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT volume = 0;
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

template<typename TriangleMesh>
typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
volume(const TriangleMesh& tmesh)
{
  return volume(tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the aspect ratio of a face of a given triangulated surface mesh.
  *
  * @tparam TriangleMesh a model of `HalfedgeGraph`
  * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
  *
  * @param f the face of which the aspect ratio is computed
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
  * @return the aspect ratio of `f`. The return type `FT` is a number type. It is
  * either deduced from the `geom_traits` \ref bgl_namedparameters "Named Parameters" if provided,
  * or the geometric traits class deduced from the point property map of `tmesh`.
  *
  */
template<typename TriangleMesh,
         typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
FT
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::FT
#endif
face_aspect_ratio(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
                  const TriangleMesh& tmesh,
                  const CGAL_PMP_NP_CLASS& np)
{
  CGAL_precondition(is_triangle(f, tmesh));

  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor           halfedge_descriptor;

  typedef typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type             Geom_traits;
  typedef typename Geom_traits::FT                                                  FT;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typename GetVertexPointMap<TriangleMesh, CGAL_PMP_NP_CLASS>::const_type
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

template<typename TriangleMesh>
typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
CGAL::vertex_point_t>::type>::Kernel::FT
face_aspect_ratio(typename boost::graph_traits<TriangleMesh>::face_descriptor f,
             const TriangleMesh& tmesh)
{
  return face_aspect_ratio(f, tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

/**
  * \ingroup measure_grp
  * computes the centroid of a volume bounded by
  * a closed triangulated surface mesh.
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
template<typename TriangleMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
#ifdef DOXYGEN_RUNNING
Point_3
#else
typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type::Point_3
#endif
centroid(const TriangleMesh& tmesh, const CGAL_PMP_NP_CLASS& np)
{
  // See: http://www2.imperial.ac.uk/~rn/centroid.pdf

  CGAL_assertion(is_triangle_mesh(tmesh));
  CGAL_assertion(is_closed(tmesh));

  using parameters::choose_parameter;
  using parameters::get_parameter;

  typedef typename GetVertexPointMap<TriangleMesh, CGAL_PMP_NP_CLASS>::const_type Vpm;
  Vpm vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_const_property_map(CGAL::vertex_point, tmesh));

  typedef typename GetGeomTraits<TriangleMesh, CGAL_PMP_NP_CLASS>::type Kernel;
  typedef typename Kernel::Point_3                                      Point_3;
  typedef typename Kernel::Vector_3                                     Vector_3;
  typedef typename Kernel::Construct_translated_point_3                 Construct_translated_point_3;
  typedef typename Kernel::Construct_vector_3                           Construct_vector_3;
  typedef typename Kernel::Construct_normal_3                           Construct_normal_3;
  typedef typename Kernel::Compute_scalar_product_3                     Scalar_product;
  typedef typename Kernel::Construct_scaled_vector_3                    Scale;
  typedef typename Kernel::Construct_sum_of_vectors_3                   Sum;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor   face_descriptor;
  typedef typename Kernel::FT FT;

  FT volume = 0;

  Vector_3 centroid(NULL_VECTOR);

  Construct_translated_point_3 point;
  Construct_vector_3 vector;
  Construct_normal_3 normal;
  Scalar_product scalar_product;
  Scale scale;
  Sum sum;

  for(face_descriptor fd : faces(tmesh))
  {
    const Point_3& p = get(vpm, target(halfedge(fd, tmesh), tmesh));
    const Point_3& q = get(vpm, target(next(halfedge(fd, tmesh), tmesh), tmesh));
    const Point_3& r = get(vpm, target(prev(halfedge(fd, tmesh), tmesh), tmesh));
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

template<typename TriangleMesh>
typename CGAL::Kernel_traits<typename property_map_value<TriangleMesh,
CGAL::vertex_point_t>::type>::Kernel::Point_3
centroid(const TriangleMesh& tmesh)
{
  return centroid(tmesh, CGAL::Polygon_mesh_processing::parameters::all_default());
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_POLYGON_MESH_PROCESSING_MEASURE_H
