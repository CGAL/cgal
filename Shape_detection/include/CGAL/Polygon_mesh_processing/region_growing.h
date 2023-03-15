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
// Author(s)     : Sébastien Loriot


#ifndef CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H
#define CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H

#include <CGAL/license/Polygon_mesh_processing/miscellaneous.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Shape_detection/Region_growing/Region_growing.h>
#include <CGAL/Shape_detection/Region_growing/Polygon_mesh.h>
#include <CGAL/Shape_detection/Region_growing/Segment_set.h>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
  \ingroup PkgPolygonMeshProcessingRef
  \brief applies a region growing algorithm to fit planes on faces of a mesh.

  See Section \ref Shape_detection_RegionGrowing for more details on the method.

  @tparam PolygonMesh a model of `FaceListGraph`
  @tparam RegionMap a model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::size_t` as value type.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param mesh the polygon mesh whose faces are used for region growing
  @param region_map a property map storing the region index of each face. Values start at `0` up to the value returned minus `1`.
         `std::size_t(-1)` is put for faces with no region assigned (can only happen if `minimum_region_size > 1`).
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  @return the number of regions detected

  \cgalNamedParamsBegin
    \cgalParamNBegin{maximum_distance}
      \cgalParamDescription{the maximum distance from a face to a plane such that it is considered part of the region of the plane}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{1}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{the maximum angle (in degrees) between the normals of the supporting planes of two adjacent faces
                            such that they are considered part of the same region}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{25 degrees}
      \cgalParamExtra{this parameter and `cosine_of_maxium_angle` are exclusive}
    \cgalParamNEnd
    \cgalParamNBegin{cosine_of_maxium_angle}
      \cgalParamDescription{The maximum angle, given as a cosine, between the normal of the supporting planes of adjacent faces
                            such that they are considered part of the same region}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{`cos(25 * PI / 180)`}
      \cgalParamExtra{this parameter and `maximum_angle` are exclusive}
    \cgalParamNEnd
    \cgalParamNBegin{minimum_region_size}
      \cgalParamDescription{the minimum number of faces such that a new region can be created from a set of faces}
      \cgalParamType{`std::size_t`}
      \cgalParamDefault{1}
    \cgalParamNEnd
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `mesh`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd
 */
template<typename PolygonMesh,
                  typename RegionMap,
                  typename NamedParameters = parameters::Default_named_parameters>
std::size_t
region_growing_of_planes_on_faces(const PolygonMesh& mesh,
                                  RegionMap region_map,
                                  const NamedParameters& np = parameters::default_values())
{
  namespace RG_PM = CGAL::Shape_detection::Polygon_mesh;

  using VPM = typename GetVertexPointMap < PolygonMesh, NamedParameters>::const_type;
  using Traits = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;

  using parameters::choose_parameter;
  using parameters::get_parameter;


  using Neighbor_query = RG_PM::One_ring_neighbor_query<PolygonMesh>;
  using Region_type = RG_PM::Least_squares_plane_fit_region<Traits, PolygonMesh, VPM>;
  using Sorting = RG_PM::Least_squares_plane_fit_sorting<Traits, PolygonMesh, Neighbor_query, VPM>;
  using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type, RegionMap>;

  Neighbor_query neighbor_query(mesh);
  Region_type region_type(mesh, np);
  Sorting sorting(mesh, neighbor_query, np);
  sorting.sort();

  std::vector<typename Region_growing::Primitive_and_region> regions;
  Region_growing region_growing(
    faces(mesh), sorting.ordered(), neighbor_query, region_type, region_map);
  region_growing.detect(CGAL::Emptyset_iterator());

  return region_growing.number_of_regions_detected();
}

/*!
  \ingroup PkgPolygonMeshProcessingRef
    \brief detects the corners on the boundary of almost planar regions by applying the region growing algorithm fitting lines on segment edges of a partition of a mesh.

    More precisely, a corner on the boundary of a region is a vertex that is either shared by at least three regions (two if it is also a vertex on the boundary of the mesh), or that is incident to two segments edges assigned to different lines.
  See Section \ref Shape_detection_RegionGrowing for more details on the method.

  @tparam PolygonMesh a model of `FaceListGraph` and `EdgeListGaph`
  @tparam RegionMap a model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%face_descriptor` as key type and `std::size_t` as value type.
  @tparam CornerIdMap a model of `WritablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor` as key type and `std::size_t` as value type.
  @tparam NamedParameters a sequence of \ref bgl_namedparameters "Named Parameters"

  @param mesh polygon mesh for region growing.
  @param region_map property map providing the region index of each face, values must be in `[0, nb_regions-1]`.
  @param corner_id_map property map storing the corner index of each vertex. Values start at `0` up to the value returned minus 1.
         `std::size_t(-1)` is put for vertices that are not corners.
  @param nb_regions the number of patches in the partition of `mesh` defined by `region_map`
  @param np an optional sequence of \ref bgl_namedparameters "Named Parameters" among the ones listed below

  @return the number of corners detected

  \cgalNamedParamsBegin
    \cgalParamNBegin{edge_is_constrained_map}
      \cgalParamDescription{a property map filled by this function such that an edge is marked as constrained
                            if it is at the interface of two different regions or on the boundary of the mesh}
      \cgalParamType{a class model of `ReadWritePropertyMap` with `boost::graph_traits<PolygonMesh>::%edge_descriptor`
                     as key type and `bool` as value type}
      \cgalParamDefault{Unused if not provided}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_distance}
      \cgalParamDescription{the maximum distance from a point to a line such that it is considered part of the region of the line}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{1}
      \cgalParamExtra{this parameter and `cosine_of_maxium_angle` are exclusive}
    \cgalParamNEnd
    \cgalParamNBegin{maximum_angle}
      \cgalParamDescription{the maximum angle in degrees between two adjacent segments
                            such that they are considered part of the same region}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{25 degrees}
    \cgalParamNEnd
    \cgalParamNBegin{cosine_of_maxium_angle}
      \cgalParamDescription{The maximum angle, given as a cosine, between two adjacent segments
                            such that they are considered part of the same region}
      \cgalParamType{`GeomTraits::FT` with `GeomTraits` being the type of the parameter `geom_traits`}
      \cgalParamDefault{`cos(25 * PI / 180)`}
      \cgalParamExtra{this parameter and `maximum_angle` are exclusive}
    \cgalParamNEnd
    \cgalParamNBegin{vertex_point_map}
      \cgalParamDescription{a property map associating points to the vertices of `mesh`}
      \cgalParamType{a class model of `ReadablePropertyMap` with `boost::graph_traits<PolygonMesh>::%vertex_descriptor`
                     as key type and `%Point_3` as value type}
      \cgalParamDefault{`boost::get(CGAL::vertex_point, mesh)`}
      \cgalParamExtra{If this parameter is omitted, an internal property map for `CGAL::vertex_point_t`
                      must be available in `PolygonMesh`.}
    \cgalParamNEnd
    \cgalParamNBegin{geom_traits}
      \cgalParamDescription{an instance of a geometric traits class}
      \cgalParamType{a class model of `Kernel`}
      \cgalParamDefault{a \cgal Kernel deduced from the point type, using `CGAL::Kernel_traits`}
      \cgalParamExtra{The geometric traits class must be compatible with the vertex point type.}
    \cgalParamNEnd
  \cgalNamedParamsEnd
 */
template <typename PolygonMesh,
          typename RegionMap,
          typename CornerIdMap,
          typename NamedParameters = parameters::Default_named_parameters>
std::size_t
detect_corners_of_regions(
  const PolygonMesh& mesh,
  RegionMap region_map,
  std::size_t nb_regions,
  CornerIdMap corner_id_map,
  const NamedParameters& np = parameters::default_values())
{
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  using Traits = typename GetGeomTraits<PolygonMesh, NamedParameters>::type;
  using Graph_traits = boost::graph_traits<PolygonMesh>;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using edge_descriptor = typename Graph_traits::edge_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;

  using Default_ecm = typename boost::template property_map<PolygonMesh, CGAL::dynamic_edge_property_t<bool> >::type;
  using Ecm = typename internal_np::Lookup_named_param_def <
                internal_np::edge_is_constrained_t,
                NamedParameters,
                Default_ecm
              > ::type;

  Default_ecm dynamic_ecm;
  if(!(is_default_parameter<NamedParameters, internal_np::edge_is_constrained_t>::value))
    dynamic_ecm = get(CGAL::dynamic_edge_property_t<bool>(), mesh);
  Ecm ecm = choose_parameter(get_parameter(np, internal_np::edge_is_constrained), dynamic_ecm);

  using Polyline_graph     = CGAL::Shape_detection::Polygon_mesh::Polyline_graph<PolygonMesh>;
  using Segment_map        = typename Polyline_graph::Segment_map;
  using Item               = typename Polyline_graph::Item;

  using Line_region  = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_region<Traits, Item, Segment_map>;
  using Line_sorting = CGAL::Shape_detection::Segment_set::Least_squares_line_fit_sorting<Traits, Item, Polyline_graph, Segment_map>;
  using RG_lines     = CGAL::Shape_detection::Region_growing<Polyline_graph, Line_region>;

  // mark as constrained edges at the interface of two regions
  for (edge_descriptor e : edges(mesh))
  {
    halfedge_descriptor h = halfedge(e, mesh);
    face_descriptor f1 = face(h, mesh);
    face_descriptor f2 = face(opposite(h, mesh), mesh);
    // an edge is constrained if it is a border edge or if the regions assigned to the two faces are different
    if (f1 == Graph_traits::null_face() || f2 == Graph_traits::null_face() || get(region_map,f1)!=get(region_map,f2))
      put(ecm, e, true);
  }

  // filter trivial edges: incident to a plane with only one face
  // such an edge cannot be removed and its vertices are corners
  std::vector<int> nb_faces_per_patch(nb_regions,0);
  for(face_descriptor f : faces(mesh))
  {
    std::size_t pid = get(region_map, f);
    nb_faces_per_patch[pid]+=1;
  }

  std::vector<edge_descriptor> filtered_edges, trivial_edges;
  for (edge_descriptor e : edges(mesh))
  {
    halfedge_descriptor h=halfedge(e,mesh);
    std::size_t r1 = is_border(h, mesh)?std::size_t(-1):get(region_map, face(h, mesh));
    h=opposite(h, mesh);
    std::size_t r2 = is_border(h, mesh)?std::size_t(-1):get(region_map, face(h, mesh));
    if ( (r1!=std::size_t(-1) && nb_faces_per_patch[r1]==1) || (r2!=std::size_t(-1) && nb_faces_per_patch[r2]==1) )
      trivial_edges.push_back(e);
    else
      filtered_edges.push_back(e);
  }

  Polyline_graph pgraph(mesh, filtered_edges, region_map);
  const auto& segment_range = pgraph.segment_range();

  Line_region line_region(np.segment_map(pgraph.segment_map()));

  Line_sorting line_sorting(
    segment_range, pgraph, CGAL::parameters::segment_map(pgraph.segment_map()));
  line_sorting.sort();

  RG_lines rg_lines(
    segment_range, pgraph, line_region);

  std::vector< std::pair<typename Line_region::Primitive, std::vector<edge_descriptor> > > subregions;
  rg_lines.detect(std::back_inserter(subregions));

#ifdef CGAL_DEBUG_DETECT_CORNERS_OF_REGIONS
  std::ofstream debug_corners("corners.xyz");
  debug_corners.precision(17);
  std::ofstream debug_edges("contraints.polylines.txt");
  debug_edges.precision(17);
#endif

  // detect vertex corner id
  std::size_t cid=0;
  for (const std::pair<typename Line_region::Primitive, std::vector<edge_descriptor>>& r : subregions)
  {
    std::vector<std::size_t> vertex_count(num_vertices(mesh), 0);
    std::vector<vertex_descriptor> line_vertices;
    auto register_vertex = [&vertex_count, &line_vertices]
                           (vertex_descriptor v)
    {
      if (vertex_count[v]==0)
        line_vertices.push_back(v);
      vertex_count[v]+=1;
    };
    for (edge_descriptor e : r.second)
    {
      put(ecm, e, true);
      register_vertex(source(e, mesh));
      register_vertex(target(e, mesh));
#ifdef CGAL_DEBUG_DETECT_CORNERS_OF_REGIONS
      debug_edges << "2 " << mesh.point(source(e, mesh)) << " " << mesh.point(target(e, mesh)) << "\n";
#endif
    }

    for (vertex_descriptor v : line_vertices)
      if (vertex_count[v]==1)
      {
#ifdef CGAL_DEBUG_DETECT_CORNERS_OF_REGIONS
        debug_corners << mesh.point(v) << "\n";
#endif
        if (get(corner_id_map, v) == std::size_t(-1))
          put(corner_id_map, v, cid++);
      }
  }

  // process trivial edges (could be done before if needed)
  for(edge_descriptor e : trivial_edges)
  {
#ifdef CGAL_DEBUG_DETECT_CORNERS_OF_REGIONS
    debug_edges << "2 " << mesh.point(source(e, mesh)) << " " << mesh.point(target(e, mesh)) << "\n";
#endif
    put(ecm, e, true);
    if (get(corner_id_map, source(e, mesh))==std::size_t(-1))
      put(corner_id_map, source(e, mesh), cid++);
    if (get(corner_id_map, target(e, mesh))==std::size_t(-1))
      put(corner_id_map, target(e, mesh), cid++);
  }

  return cid;
}

} } // end of CGAL::Polygon_mesh_processing namespace

#endif //CGAL_POLYGON_MESH_PROCESSING_REGION_GROWING_H
