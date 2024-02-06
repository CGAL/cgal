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

#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/boost/graph/copy_face_graph.h>

#include <boost/graph/graph_traits.hpp>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/boost/graph/named_params_helper.h>

#include <vector>

namespace CGAL {

  namespace Polygon_mesh_processing {

    /*!
    * \ingroup PkgPolygonMeshProcessingRef
    *
    * computes a bounding box of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `VertexListGraph`
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
    *                           `Construct_bbox_3` must provide the functor `Bbox_3 operator()(Point_3)`
    *                           where `%Point_3` is the value type of the vertex point map.}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @see `vertex_bbox()`
    * @see `edge_bbox()`
    * @see `face_bbox()`
    */
    template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
    CGAL::Bbox_3 bbox(const PolygonMesh& pmesh,
                      const NamedParameters& np = parameters::default_values())
    {
      using CGAL::parameters::choose_parameter;
      using CGAL::parameters::get_parameter;

      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
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
    *
    * computes a bounding box of the vertex of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `Graph`
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
    * @see `edge_bbox()`
    * @see `face_bbox()`
    * @see `bbox()`
    */
    template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
    CGAL::Bbox_3 vertex_bbox(typename boost::graph_traits<PolygonMesh>::vertex_descriptor vd,
                             const PolygonMesh& pmesh,
                             const NamedParameters& np = parameters::default_values())
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
    *
    * computes a bounding box of an edge of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `Graph`
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
    * @see `vertex_bbox()`
    * @see `face_bbox()`
    * @see `bbox()`
    */
    template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
    CGAL::Bbox_3 edge_bbox(typename boost::graph_traits<PolygonMesh>::edge_descriptor ed,
                           const PolygonMesh& pmesh,
                           const NamedParameters& np = parameters::default_values())
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      CGAL_precondition(is_valid_edge_descriptor(ed, pmesh));

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
    *
    * computes a bounding box of a face of a polygon mesh.
    *
    * @tparam PolygonMesh a model of `Graph`
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
    * @see `vertex_bbox()`
    * @see `edge_bbox()`
    * @see `bbox()`
    */
    template<typename PolygonMesh, typename NamedParameters = parameters::Default_named_parameters>
    CGAL::Bbox_3 face_bbox(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                           const PolygonMesh& pmesh,
                           const NamedParameters& np = parameters::default_values())
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      CGAL_precondition(is_valid_face_descriptor(fd, pmesh));

      typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type
        vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(CGAL::vertex_point, pmesh));

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_bbox_3 get_bbox = gt.construct_bbox_3_object();

      typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;

      CGAL::Bbox_3 bb;
      for(halfedge_descriptor h : halfedges_around_face(halfedge(fd, pmesh), pmesh))
      {
        bb += get_bbox( get(vpm, target(h, pmesh)) );
      }
      return bb;
    }

    /*!
    * adds a triangulated bounding box to a polygon mesh
    * @todo add extended bbox factor
    * @todo add triangulate or not as NP
    */
    template<typename PolygonMesh,
             typename NamedParameters = parameters::Default_named_parameters>
    void add_bbox(PolygonMesh& pmesh,
                  const NamedParameters& np = parameters::default_values())
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GT;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typedef typename GT::Point_3 Point_3;

      const double factor = 1.1; //todo : add as NP
      CGAL::Bbox_3 bb = bbox(pmesh);
      if (factor != 1.0)
      {
        const double dx = bb.xmax() - bb.xmin();
        const double dy = bb.ymax() - bb.ymin();
        const double dz = bb.zmax() - bb.zmin();
        const Point_3 center( bb.xmin() + 0.5 * dx,
                              bb.ymin() + 0.5 * dy,
                              bb.zmin() + 0.5 * dz );
        bb = Bbox_3( center.x() - factor * 0.5 * dx,
                     center.y() - factor * 0.5 * dy,
                     center.z() - factor * 0.5 * dz,
                     center.x() + factor * 0.5 * dx,
                     center.y() + factor * 0.5 * dy,
                     center.z() + factor * 0.5 * dz );
      }

      const bool triangulate = true; //todo : add as NP

      const typename GT::Iso_cuboid_3 bbox(bb);
      std::vector<Point_3> bb_points;
      for (int i = 0; i < 8; ++i)
        bb_points.push_back(bbox[i]);

      std::vector<std::vector<std::size_t>> faces;
      if (!triangulate)
      {
        faces.push_back({0, 1, 2, 3});//bottom
        faces.push_back({4, 5, 6, 7});//top
        faces.push_back({0, 1, 6, 5});//front
        faces.push_back({2, 3, 4, 7});//back
        faces.push_back({1, 2, 7, 6});//right
        faces.push_back({0, 3, 4, 5});//left
      }
      else
      {
        faces.push_back({0, 1, 2});//bottom
        faces.push_back({0, 2, 3});
        faces.push_back({4, 5, 6});//top
        faces.push_back({4, 6, 7});
        faces.push_back({0, 1, 5});//front
        faces.push_back({1, 5, 6});
        faces.push_back({2, 3, 4});//back
        faces.push_back({2, 4, 7});
        faces.push_back({1, 2, 6});//right
        faces.push_back({2, 6, 7});
        faces.push_back({0, 3, 4});//left
        faces.push_back({0, 4, 5});
      }

      PolygonMesh bbox_pmesh;
      orient_polygon_soup(bb_points, faces);
      polygon_soup_to_polygon_mesh(bb_points, faces, bbox_pmesh);

      CGAL::copy_face_graph(bbox_pmesh, pmesh,
                            parameters::default_values(),
                            np);
    }
  }
}

#endif //CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H
