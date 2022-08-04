// Copyright (c) 2022 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H
#define CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
  \ingroup PkgPolygonMeshProcessingRef
  applies the region growing algorithm to fit planes on faces of `mesh`.
  See Section \ref Shape_detection_RegionGrowing for more details on the method.

  @tparam PolygonMesh
    a model of `FaceListGraph`
  @tparam RegionMap a model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_iterator` as key type and `std::size_t` as value_type.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param mesh polygon mesh for region growing.
  @param region_map property map storing the region index of each face. Values start at `0` up to the value returned minus 1.
         `std::size_t(-1)` is put for faces with no region assigned (can only append if minimum_region_size > 1).
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  @return the number of regions detected

  \cgalNamedParamsBegin
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `pmesh`}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, pmesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_distance}
      \cgalParamDescription{the maximum distance from a point to a line}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{1}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{the maximum angle in degrees between
      the normal of a point and the normal of a line}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{cosine_value}
      \cgalParamDescription{the cos value computed as `cos(maximum_angle * PI / 180)`,
      this parameter can be used instead of the `maximum_angle`}
      \cgalParamType{`GeomTraits::FT`}
      \cgalParamDefault{`cos(25 * PI / 180)`}
    \cgalParamNEnd
    \cgalParamNBegin{minimum_region_size}
      \cgalParamDescription{the minimum number of points a region must have}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{1}
    \cgalParamNEnd
  \cgalNamedParamsEnd
 */
template<
typename PolygonMesh,
typename RegionMap,
typename CGAL_NP_TEMPLATE_PARAMETERS>
std::size_t
region_growing_of_planes_on_faces(
  const PolygonMesh& mesh,
  RegionMap region_map,
  const CGAL_NP_CLASS& np = parameters::default_values())
{
  typedef typename GetVertexPointMap < PolygonMesh, CGAL_NP_CLASS>::const_type VPM;
  typedef typename GetGeomTraits<PolygonMesh, CGAL_NP_CLASS>::type  Traits;

  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::face_descriptor face_descriptor;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  namespace RG_PM = CGAL::Shape_detection::Polygon_mesh;

  using Neighbor_query = RG_PM::One_ring_neighbor_query<PolygonMesh>;
  using Region_type = RG_PM::Least_squares_plane_fit_region<Traits, PolygonMesh, VPM>;
  using Sorting = RG_PM::Least_squares_plane_fit_sorting<Traits, PolygonMesh, Neighbor_query, VPM>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;

  Neighbor_query neighbor_query(mesh);
  Region_type region_type(mesh, np);
  Sorting sorting(mesh, neighbor_query, np);
  sorting.sort();

  std::vector<typename Region_growing::Primitive_and_region> regions;
  Region_growing region_growing(
    faces(mesh), sorting.ordered(), neighbor_query, region_type);
  region_growing.detect(CGAL::Emptyset_iterator());

  for (face_descriptor f : faces(mesh))
    put(region_map, f, get(region_growing.region_map(), f));

  return region_growing.number_of_regions_detected();
}

} } // end of CGAL::Polygon_mesh_processing namespace

#endif //CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H
