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
// Author(s)     : Jane Tournois

#ifndef CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H
#define CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>


#include <CGAL/Bbox_3.h>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>


#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

#ifdef DOXYGEN_RUNNING
#define CGAL_PMP_NP_TEMPLATE_PARAMETERS NamedParameters
#define CGAL_PMP_NP_CLASS NamedParameters
#endif

namespace CGAL {

  namespace Polygon_mesh_processing {

    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    *  computes a bounding box of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `HalfedgeListGraph`
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param pmesh a polygon mesh
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
    *                     must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{geom_traits}
    *     \cgalParamDescription{an instance of a geometric traits class providing the functor `Construct_bbox_3`
    *                           and the function `Construct_bbox_3 construct_bbox_3_object()`.
    *                           `Construct_bbox_3` must provide `Bbox_3 operator()(Point_3)`
    *                           where `%Point_3` is the value type of the vertex point map.}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @return a bounding box of `pmesh`
    */
    template<typename PolygonMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
    CGAL::Bbox_3 bbox(const PolygonMesh& pmesh,
                      const CGAL_PMP_NP_CLASS& np)
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      typename GetVertexPointMap<PolygonMesh, CGAL_PMP_NP_CLASS>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, CGAL_PMP_NP_CLASS>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_bbox_3 get_bbox = gt.construct_bbox_3_object();

      typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

      CGAL::Bbox_3 bb;
      for(vertex_descriptor v : vertices(pmesh))
      {
        bb += get_bbox( get(vpm, v) );
      }
      return bb;
    }

    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    *  computes a bounding box of a vertex of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `HalfedgeGraph`
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param vd a descriptor of a vertex in `pmesh`
    * @param pmesh a polygon mesh
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
    *                     must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{geom_traits}
    *     \cgalParamDescription{an instance of a geometric traits class providing the functor `Construct_bbox_3`
    *                           and the function `Construct_bbox_3 construct_bbox_3_object()`.
    *                           `Construct_bbox_3` must provide `Bbox_3 operator()(Point_3)`
    *                           where `%Point_3` is the value type of the vertex point map.}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @return a bounding box of `pmesh`
    */
    template<typename PolygonMesh, typename NamedParameters>
    CGAL::Bbox_3 vertex_bbox(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                             const PolygonMesh& pmesh,
                             const NamedParameters& np)
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;
      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_bbox_3 get_bbox = gt.construct_bbox_3_object();

      return get_bbox( get(vpm, vd) );
    }

    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    *  computes a bounding box of an edge of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `HalfedgeGraph`
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param ed a descriptor of an edge in `pmesh`
    * @param pmesh a polygon mesh
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
    *                     must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{geom_traits}
    *     \cgalParamDescription{an instance of a geometric traits class providing the functor `Construct_bbox_3`
    *                           and the function `Construct_bbox_3 construct_bbox_3_object()`.
    *                           `Construct_bbox_3` must provide `Bbox_3 operator()(Point_3)`
    *                           where `%Point_3` is the value type of the vertex point map.}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @return a bounding box of `pmesh`
    */
    template<typename PolygonMesh, typename NamedParameters>
    CGAL::Bbox_3 edge_bbox(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed,
                           const PolygonMesh& pmesh,
                           const NamedParameters& np)
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_bbox_3 get_bbox = gt.construct_bbox_3_object();

      return get_bbox( get(vpm, source(ed, pmesh)) ) +
             get_bbox( get(vpm, target(ed, pmesh)) );
    }

    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    *  computes a bounding box of a face of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `HalfedgeGraph`
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param fd a descriptor of a face in `pmesh`
    * @param pmesh a polygon mesh
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
    *                     must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *
    *   \cgalParamNBegin{geom_traits}
    *     \cgalParamDescription{an instance of a geometric traits class providing the functor `Construct_bbox_3`
    *                           and the function `Construct_bbox_3 construct_bbox_3_object()`.
    *                           `Construct_bbox_3` must provide `Bbox_3 operator()(Point_3)`
    *                           where `%Point_3` is the value type of the vertex point map.}
    *     \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @return a bounding box of `pmesh`
    */
    template<typename PolygonMesh, typename NamedParameters>
    CGAL::Bbox_3 face_bbox(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                           const PolygonMesh& pmesh,
                           const NamedParameters& np)
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_bbox_3 get_bbox = gt.construct_bbox_3_object();

      typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

      CGAL::Bbox_3 bb;
      for(halfedge_descriptor h :
                    halfedges_around_face(halfedge(fd, pmesh), pmesh))
      {
        bb += get_bbox( get(vpm, target(h, pmesh)) );
      }
      return bb;
    }

    template<typename PolygonMesh>
    CGAL::Bbox_3 vertex_bbox(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                               const PolygonMesh& pmesh)
    {
      return vertex_bbox(vd, pmesh,
        CGAL::Polygon_mesh_processing::parameters::all_default());
    }
    template<typename PolygonMesh>
    CGAL::Bbox_3 edge_bbox(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed,
                             const PolygonMesh& pmesh)
    {
      return edge_bbox(ed, pmesh,
        CGAL::Polygon_mesh_processing::parameters::all_default());
    }
    template<typename PolygonMesh>
    CGAL::Bbox_3 face_bbox(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                             const PolygonMesh& pmesh)
    {
      return face_bbox(fd, pmesh,
        CGAL::Polygon_mesh_processing::parameters::all_default());
    }

    template<typename PolygonMesh>
    CGAL::Bbox_3 bbox(const PolygonMesh& pmesh)
    {
      return bbox(pmesh,
        CGAL::Polygon_mesh_processing::parameters::all_default());
    }

    // deprecated function
    #ifndef CGAL_NO_DEPRECATED_CODE
    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    * \deprecated This function is deprecated since \cgal 4.10, `CGAL::Polygon_mesh_processing::bbox()` should be used instead.
    */
    template<typename PolygonMesh, typename CGAL_PMP_NP_TEMPLATE_PARAMETERS>
    CGAL_DEPRECATED CGAL::Bbox_3 bbox_3(const PolygonMesh& pmesh,
                        const CGAL_PMP_NP_CLASS& np)
    {
      return bbox(pmesh, np);
    }

    template<typename PolygonMesh>
    CGAL_DEPRECATED CGAL::Bbox_3 bbox_3(const PolygonMesh& pmesh)
    {
      return bbox(pmesh);
    }
    #endif // CGAL_NO_DEPRECATED_CODE
  }
}

#endif //CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H

