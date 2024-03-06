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

#include <CGAL/boost/graph/generators.h>

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
    *
    *   \cgalParamNBegin{bbox_scaling}
    *     \cgalParamDescription{a double used to scale the bounding box.
    *       The default value is 1 and the bounding box is the smallest possible
    *       axis-aligned bounding box.}
    *     \cgalParamDefault{1.}
    *     \cgalParamPrecondition{`bbox_scaling > 0`}
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

      const double factor = choose_parameter(get_parameter(np, internal_np::bbox_scaling), 1.);
      CGAL_precondition(factor > 0);

      typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

      CGAL::Bbox_3 bb;
      for(vertex_descriptor v : vertices(pmesh))
      {
        bb += get_bbox( get(vpm, v) );
      }
      bb.scale(factor);

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
    * \ingroup PkgPolygonMeshProcessingRef
    *
    * adds an axis-aligned bounding box to a polygon mesh.
    *
    * @tparam PolygonMesh a model of `MutableFaceGraph`
    * @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"
    *
    * @param pmesh a polygon mesh
    * @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below
    *
    * \cgalNamedParamsBegin
    *   \cgalParamNBegin{bbox_scaling}
    *     \cgalParamDescription{a double used to scale the bounding box.
    *       The default value is 1 and the bounding box is the smallest possible
    *       axis-aligned bounding box.}
    *     \cgalParamDefault{1.}
    *     \cgalParamPrecondition{`bbox_scaling > 0`}
    *   \cgalParamNEnd
    *   \cgalParamNBegin{do_not_triangulate_faces}
    *     \cgalParamDescription{a Boolean used to specify whether the bounding box's faces
    *       should be triangulated or not.
    *       The default value is `true`, and faces are not triangulated.}
    *     \cgalParamDefault{true}
    *   \cgalParamNEnd
    *   \cgalParamNBegin{vertex_point_map}
    *     \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
    *     \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
    *                    as key type and `%Point_3` as value type}
    *     \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
    *     \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
    *                     must be available in `PolygonMesh`.}
    *   \cgalParamNEnd
    *   \cgalParamNBegin{geom_traits}
    *     \cgalParamDescription{an instance of a geometric traits class model of `Kernel`.}
    *   \cgalParamNEnd
    * \cgalNamedParamsEnd
    *
    * @see `bbox()`
    */
    template<typename PolygonMesh,
             typename NamedParameters = parameters::Default_named_parameters>
    void add_bbox(PolygonMesh& pmesh,
                  const NamedParameters& np = parameters::default_values())
    {
      using parameters::choose_parameter;
      using parameters::get_parameter;

      using GT = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
      using P = typename GT::Point_3;
      GT gt = choose_parameter<GT>(get_parameter(np, internal_np::geom_traits));
      typename GT::Construct_iso_cuboid_3
        iso_cuboid = gt.construct_iso_cuboid_3_object();

      const CGAL::Bbox_3 bb = bbox(pmesh, np);
      CGAL::make_hexahedron(iso_cuboid(P(bb.xmin(), bb.ymin(), bb.zmin()),
                                       P(bb.xmax(), bb.ymax(), bb.zmax())),
                            pmesh, np);
    }
  }
}

#endif //CGAL_POLYGON_MESH_PROCESSING_BOUNDING_BOX_H
