// Copyright (c) 2024-2025 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H
#define CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H

#include <CGAL/license/Straight_skeleton_3.h>

#include <CGAL/Straight_skeleton_3/Configuration.h>
#include <CGAL/Straight_skeleton_3/internal/HDS/Polyhedron.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/HDS_utils.h>
#include <CGAL/Straight_skeleton_3/internal/algorithm/Polyhedron_transformation.h>

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/mark_domain_in_triangulation.h>
#include <CGAL/property_map.h>
#if 0
# include <CGAL/Polygon_mesh_processing/region_growing.h>
#endif
# include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/unordered_flat_map.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <unordered_map>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace IO {

template <typename GeomTraits>
class FaceGraphIO
{
  using FT = typename GeomTraits::FT;
  using Point_3 = typename GeomTraits::Point_3;
  using Vector_3 = typename GeomTraits::Vector_3;
  using Plane_3 = typename GeomTraits::Plane_3;

private:
  using Polyhedron = internal::HDS::Polyhedron<GeomTraits>;
  using PolyhedronSPtr = typename Polyhedron::PolyhedronSPtr;

  using Vertex = typename Polyhedron::Vertex;
  using VertexSPtr = typename Polyhedron::VertexSPtr;
  using EdgeWPtr = typename Polyhedron::EdgeWPtr;
  using EdgeSPtr = typename Polyhedron::EdgeSPtr;
  using Facet = typename Polyhedron::Facet;
  using FacetSPtr = typename Polyhedron::FacetSPtr;

  using Skeleton_facet_data = typename Polyhedron::Skeleton_facet_data;
  using SkelFacetDataSPtr = typename Polyhedron::SkelFacetDataSPtr;

private:
  using Transformation = internal::algorithm::Polyhedron_transformation<GeomTraits>;
  using Hds_utils = internal::algorithm::Hds_utils<GeomTraits>;

public:
  template <typename TriangleMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr load(const TriangleMesh& tmesh,
                             std::map<typename boost::graph_traits<TriangleMesh>::edge_descriptor, EdgeWPtr>& e2e,
                             const NamedParameters& np = CGAL::parameters::default_values())
  {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using vertex_descriptor = typename boost::graph_traits<TriangleMesh>::vertex_descriptor;
    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    using VPM = typename GetVertexPointMap<TriangleMesh, NamedParameters>::const_type;
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_const_property_map(vertex_point, tmesh));

    const bool outward_offsetting = choose_parameter(get_parameter(np, internal_np::outward_offsetting), false);

    CGAL_warning(CGAL::is_triangle_mesh(tmesh));

    PolyhedronSPtr result = Polyhedron::create();

    unsigned int vertex_id_new = 0;

    for (vertex_descriptor vi : vertices(tmesh)) {
      ++vertex_id_new;
      decltype(auto) point = get(vpm, vi);
      VertexSPtr vertex = Vertex::create(point);
      vertex->set_id(vertex_id_new);
      result->add_vertex(vertex);
    }

    std::vector<VertexSPtr> vertices(result->vertices().begin(), result->vertices().end());

    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, FT>(1));

    int facet_id_new = -1;
    for (face_descriptor fi : faces(tmesh)) {
      ++facet_id_new;

      unsigned int num_vertices = degree(fi, tmesh);
      CGAL_SS3_IO_TRACE_V(16, "new: F" << facet_id_new << " with " << num_vertices << " vertices");
      CGAL_assertion(num_vertices > 2);

      std::vector<VertexSPtr> poly_vertices(num_vertices);
      for (unsigned int i = 0; i < num_vertices; ++i) {
        poly_vertices[i] = VertexSPtr();
      }

      unsigned int pos = 0;
      for (halfedge_descriptor h : halfedges_around_face(halfedge(fi, tmesh), tmesh)) {
        unsigned int vertex_id = source(h, tmesh);
        if (vertex_id < vertices.size()) {
          poly_vertices[pos++] = vertices[vertex_id];
          CGAL_SS3_IO_TRACE_V(32, "  V" << vertices[vertex_id]->id() << "; "
                                        << vertices[vertex_id]->point());
        } else {
          std::stringstream whatstream;
          whatstream << "Vertex with id=" << vertex_id << " does not exist.";
          throw std::runtime_error(whatstream.str());
        }
      }

      if (outward_offsetting)
        std::reverse(poly_vertices.begin(), poly_vertices.end());


      FacetSPtr facet = Facet::create(poly_vertices);
      facet->set_id(facet_id_new);

      // Correspondence between the edges of the input mesh and the new edges
      // in the polyhedron
      // poly_vertices is filled starting at source() of the first edge
      // Facet::create() creates the i-th edge between vertices[i] and vertices[i+1]
      halfedge_descriptor h = halfedge(fi, tmesh);
      for (EdgeSPtr e : facet->edges()) {
        e2e[edge(h, tmesh)] = e;
        h = next(h, tmesh);
      }

      Vector_3 n = Polygon_mesh_processing::compute_face_normal(fi, tmesh);
      if (outward_offsetting)
        n = -n;

      Plane_3 pl (poly_vertices[0]->point(), n);
      facet->set_plane(pl);
      result->add_facet(facet);

      const FT weight = get(weight_pmap, fi);
      CGAL_assertion(weight >= 0);

      SkelFacetDataSPtr data = Skeleton_facet_data::create(facet);
      data->set_speed(weight);
      data->set_input_face_id(facet_id_new);
    }

    for (const EdgeSPtr& edge : result->edges()) {
      if (!(edge->get_facet_L() && edge->get_facet_R())) {
        CGAL_SS3_IO_TRACE_V(1, "Warning: Polyhedron has no closed boundary.");
        CGAL_SS3_IO_TRACE_V(1, edge->to_string());
      }
    }

    Transformation::remove_vertices_deg_lt3(result);

    CGAL_SS3_DEBUG_SPTR(result);
    CGAL_postcondition(result->is_consistent());

    return result;
  }

  template <typename TriangleMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr load(const TriangleMesh& tmesh,
                             const NamedParameters& np = CGAL::parameters::default_values())
  {
    using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
    std::map<edge_descriptor, EdgeWPtr> unused_e2e;
    return load(tmesh, unused_e2e, np);
  }

  template <typename TriangleMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr convert(const TriangleMesh& tmesh,
                                const NamedParameters& np = CGAL::parameters::default_values())
  {
    CGAL_SS3_TRANSF_TRACE("Converting mesh...");

    ConfigurationSPtr config = Configuration::get_instance();
    bool merge_faces = config->get_Boolean("Preprocessing", "merge_coplanar_faces");

    PolyhedronSPtr polyhedron = IO::FaceGraphIO<GeomTraits>::load(tmesh, np);
    if (!merge_faces)
      return polyhedron;

#if 1

#if 0
    // this actually makes things worse because we want to merge almost coplanar facets afterwards,
    // and truncating precision will increase the variance of the normals of facets almost living
    // on the same non-cardinal plane, resulting in fewer merges.
    Transformation::truncate_precision(polyhedron);
#endif

    Transformation::merge_coplanar_facets(polyhedron);
#else
    namespace PMP = CGAL::Polygon_mesh_processing;

    CGAL::Bbox_3 bbox = PMP::bbox(tmesh);
    const FT diag_length = CGAL::approximate_sqrt(square(bbox.xmax() - bbox.xmin()) +
                                                  square(bbox.ymax() - bbox.ymin()) +
                                                  square(bbox.zmax() - bbox.zmin()));

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(tmesh));
    boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected

    const FT cos_of_max_angle = 0.98;
    const FT max_distance = 0.0001 * diag_length;

    // detect planar regions in the mesh
    // @todo growing should:
    // - use the .ini value of 'coplanarity_epsilon'
    // - stop if it merges faces with different weights
    // - give an error for adjacent coplanar faces that have different weights
    std::size_t nb_regions =
        PMP::region_growing_of_planes_on_faces(tmesh,
                                               CGAL::make_random_access_property_map(region_ids),
                                               CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle)
                                                               .region_primitive_map(plane_map)
                                                               .maximum_distance(max_distance));

    static int region_dump_id = -1;
    utils::save_colored_mesh(tmesh, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");

    // detect corner vertices on the boundary of planar regions
    std::vector<std::size_t> corner_ids(num_vertices(tmesh), -1); // corner status of vertices
    std::vector<bool> ecm(num_edges(tmesh), false); // mark edges at the boundary of regions

    std::size_t nb_corners =
        PMP::detect_corners_of_regions(tmesh,
                                      CGAL::make_random_access_property_map(region_ids),
                                      nb_regions,
                                      CGAL::make_random_access_property_map(corner_ids),
                                      CGAL::parameters::cosine_of_maximum_angle(cos_of_max_angle).
                                                        maximum_distance(max_distance).
                                                        edge_is_constrained_map(CGAL::make_random_access_property_map(ecm)));

    CGAL_SS3_TRANSF_TRACE_CODE(for (face_descriptor f : faces(tmesh)))
    CGAL_SS3_TRANSF_TRACE("facet " << f << " is in region " << region_ids[f]);

    // the almost-coplanar merge is performed after the conversion to the Polyhedron
    // data structure because we want to be able to create faces that have holes,
    // which the CGAL::Surface_mesh class does not support
    std::map<edge_descriptor, EdgeWPtr> e2e;
    PolyhedronSPtr polyhedron = db::_3d::FaceGraphIO::load(sm, e2e, np);

    // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
    for (edge_descriptor e: edges(tmesh)) {
        if (ecm[e]) {
          continue;
        }

        EdgeSPtr edge = e2e[e].lock();
        if (!edge) {
          continue;
        }

        CGAL_SS3_TRANSF_TRACE("Merging facets " << edge->get_facet_L()->id() << " and " << edge->get_facet_R()->id());
        CGAL_assertion(sm.point(source(e, sm)) == edge->source()->point());
        CGAL_assertion(sm.point(target(e, sm)) == edge->target()->point());

        // @fixme it seems like intermediate states are somewhat unsound during edge merging
        merge_facets(edge, polyhedron);
    }

    polyhedron->initialize_all_IDs();

    sanitize(polyhedron);
#endif

    CGAL_SS3_TRANSF_TRACE("Converted, " << polyhedron->facets().size() << " facets");

    return polyhedron;
  }

public:
  template <typename PolygonMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static bool save(const PolyhedronSPtr& polyhedron,
                   PolygonMesh& pmesh,
                   const NamedParameters& np = CGAL::parameters::default_values())
  {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    using vertex_descriptor = typename boost::graph_traits<PolygonMesh>::vertex_descriptor;
    using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

    using Itag = CGAL::Exact_intersections_tag;
    using PK = CGAL::Projection_traits_3<GeomTraits>;
    using PVbb = CGAL::Triangulation_vertex_base_with_info_2<VertexSPtr, PK>;
    using PVb = CGAL::Triangulation_vertex_base_2<PK, PVbb>;
    using PFb = CGAL::Constrained_triangulation_face_base_2<PK>;
    using PTDS = CGAL::Triangulation_data_structure_2<PVb, PFb>;
    using PCDT = CGAL::Constrained_Delaunay_triangulation_2<PK, PTDS, Itag>;
    using PCDT_VH = typename PCDT::Vertex_handle;
    using PCDT_FH = typename PCDT::Face_handle;

    using VPM = typename GetVertexPointMap<PolygonMesh, NamedParameters>::type;
    VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                               get_property_map(vertex_point, pmesh));

    CGAL_SS3_IO_TRACE_V(8, "Save polyhedron with " << polyhedron->vertices().size() << " vertices and "
                                                   << polyhedron->facets().size() << " facets");

    // @todo do not systematically triangulate, but use this NP and if it is false,
    // only triangulate what is not representable otherwise (see code in PMP::remesh_planar_faces)
    bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

    // Vertices
    CGAL::unordered_flat_map<VertexSPtr, vertex_descriptor> v_map;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      vertex_descriptor vi = add_vertex(pmesh);
      put(vpm, vi, vertex->point());
      v_map[vertex] = vi;
    }

    // Write facets
    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, FT>(1));

    auto f2i = choose_parameter(get_parameter(np, internal_np::face_to_face_map),
                                boost::make_assoc_property_map(std::map<face_descriptor, std::size_t>()));

    for (const FacetSPtr& facet : polyhedron->facets()) {
      CGAL_assertion(facet->id() != -1);
      CGAL_assertion(facet->edges().size() >= 3);

      FT speed = Hds_utils::get_speed(facet);
      std::size_t input_face_id = Hds_utils::get_input_face_id(facet);
      CGAL_SS3_IO_TRACE_V(16, "Saving facet " << facet->id() << " with speed " << speed << " and input face id " << input_face_id);

      Vector_3 n = facet->get_plane().orthogonal_vector();
      CGAL_assertion(n != CGAL::NULL_VECTOR);

      PK traits(n);
      PCDT pcdt(traits);

      std::map<VertexSPtr, PCDT_VH> face_vhs;

      for (const VertexSPtr& vertex : facet->vertices()) {
        auto res = face_vhs.emplace(vertex, PCDT_VH());
        if (res.second) { // first time seeing this point
          PCDT_VH vh = pcdt.insert(vertex->point());
          res.first->second = vh;
          vh->info() = vertex;
        }
      }

      unsigned int ne = 0;
      for (const EdgeSPtr& edge : facet->edges()) {
        VertexSPtr v0 = edge->src(facet);
        VertexSPtr v1 = edge->tgt(facet);

        if(v0->point() == v1->point())
        {
          CGAL_SS3_IO_TRACE("Warning: degenerate edge at " << v0->point());
          break;
        }
        else
        {
          PCDT_VH vh0 = face_vhs.at(v0);
          PCDT_VH vh1 = face_vhs.at(v1);

          try
          {
            pcdt.insert_constraint(vh0, vh1);
          }
          catch(const typename PCDT::Intersection_of_constraints_exception&)
          {
            CGAL_SS3_IO_TRACE("Error: Intersection of constraints");
            CGAL_SS3_IO_TRACE("While inserting " << v0->point() << " || " << v1->point());
            CGAL_SS3_IO_TRACE(facet->to_string());
            CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
            return false;
          }
          ++ne;
        }
      }

      if(ne < 3) // degenerate facet
      {
        CGAL_SS3_IO_TRACE("Warning: skipping degenerate facet");
        continue;
      }

      std::unordered_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(pcdt, in_domain);

      for (auto fh : pcdt.finite_face_handles()) {
        if(!get(in_domain, fh))
          continue;

        auto vr = CGAL::make_array(v_map[fh->vertex(0)->info()],
                                   v_map[fh->vertex(1)->info()],
                                   v_map[fh->vertex(2)->info()]);
        face_descriptor sm_f = CGAL::Euler::add_face(vr, pmesh);
        if (sm_f == boost::graph_traits<PolygonMesh>::null_face()) {
          CGAL_SS3_IO_TRACE("Error: failed to add face to surface mesh (2)");
          CGAL_SS3_IO_TRACE("Face:\n" << fh->vertex(0)->point() << " "
                                      << fh->vertex(1)->point() << " "
                                      << fh->vertex(2)->point());
          CGAL::IO::write_polygon_mesh("results/failed.off", pmesh, CGAL::parameters::stream_precision(17));

          CGAL_assertion_code(bool cannot_add = CGAL::Euler::can_add_face(vr, pmesh, true /*verbose*/);)
          CGAL_assertion(!cannot_add);

          return false;
        }

        put(weight_pmap, sm_f, speed);
        put(f2i, sm_f, input_face_id);
      }
    }

    return true;
  }

};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H */
