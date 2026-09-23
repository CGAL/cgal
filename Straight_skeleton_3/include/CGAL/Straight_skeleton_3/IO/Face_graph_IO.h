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

// #define CGAL_SS3_DETECT_COPLANARITIES_WITH_NORMAL_CHANGE

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
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>
#ifndef CGAL_SS3_DETECT_COPLANARITIES_WITH_NORMAL_CHANGE
# include <CGAL/Polygon_mesh_processing/bbox.h>
# include <CGAL/Polygon_mesh_processing/region_growing.h>
#endif
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/unordered_flat_map.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <vector>
#include <unordered_map>

namespace CGAL {
namespace Straight_skeletons_3 {
namespace IO {
namespace utils {

// Simple helper function to draw a mesh whose faces are colored according to the weights (speeds).
template<typename PolygonMesh, typename Values>
void save_colored_mesh(const PolygonMesh& pmesh,
                       const Values& values,
                       const std::filesystem::path& fullpath)
{
  using Color = CGAL::IO::Color;

  using face_descriptor = typename boost::graph_traits<PolygonMesh>::face_descriptor;

  using value_type = typename CGAL::cpp20::remove_cvref<decltype(values[face_descriptor()])>::type;

  std::cout << "Saving colored mesh to " << fullpath << std::endl;

  // get a unique vector of values
  std::vector<value_type> unique_values;
  for (auto f : faces(pmesh)) {
    unique_values.push_back(values[f]);
  }

  std::sort(unique_values.begin(), unique_values.end());
  unique_values.erase(std::unique(unique_values.begin(), unique_values.end()), unique_values.end());

  std::cout << "Number of unique values: " << unique_values.size() << std::endl;

  std::mt19937_64 gen(0);
  std::uniform_int_distribution<> dist(0, 255);

  std::map<value_type, CGAL::Color> colors;
  for (const auto& value : unique_values) {
    colors[value] = Color(static_cast<unsigned char>(dist(gen)),
                          static_cast<unsigned char>(dist(gen)),
                          static_cast<unsigned char>(dist(gen)));
    std::cout << " value " << value << " has color " << colors[value] << std::endl;
  }

  auto& nc_pmesh = const_cast<PolygonMesh&>(pmesh);
  auto face_color = nc_pmesh.template add_property_map<face_descriptor, Color>("f:color").first;

  for (auto f : faces(pmesh)) {
    // std::cout << "facet " << f << " with value " << values[f] << " gets color " << colors[values[f]] << std::endl;
    put(face_color, f, colors[values[f]]);
  }

  std::ofstream out(fullpath);
  out.precision(17);
  CGAL::IO::write_PLY(out, pmesh, CGAL::parameters::face_color_map(face_color));
}

} // namepace utils

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
  using FacetWPtr = typename Polyhedron::FacetWPtr;
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
                             CGAL::unordered_flat_map<typename boost::graph_traits<TriangleMesh>::edge_descriptor, EdgeWPtr>& e2e,
                             CGAL::unordered_flat_map<typename boost::graph_traits<TriangleMesh>::face_descriptor, FacetWPtr>& f2f,
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

    for (vertex_descriptor vd : vertices(tmesh)) {
      ++vertex_id_new;
      decltype(auto) point = get(vpm, vd);
      VertexSPtr vertex = Vertex::create(point);
      vertex->set_id(vertex_id_new);
      result->add_vertex(vertex);
    }

    std::vector<VertexSPtr> vertices(result->vertices().begin(), result->vertices().end());

    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, FT>(1));

    int facet_id_new = -1;
    for (face_descriptor fd : faces(tmesh)) {
      ++facet_id_new;

      unsigned int num_vertices = degree(fd, tmesh);
      CGAL_SS3_IO_TRACE_V(16, "new: F" << facet_id_new << " with " << num_vertices << " vertices");
      CGAL_assertion(num_vertices > 2);

      std::vector<VertexSPtr> poly_vertices(num_vertices);
      for (unsigned int i = 0; i < num_vertices; ++i) {
        poly_vertices[i] = VertexSPtr();
      }

      unsigned int pos = 0;
      for (halfedge_descriptor h : halfedges_around_face(halfedge(fd, tmesh), tmesh)) {
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
      f2f[fd] = facet;

      // let's not assume anything on the edge order...
      std::map<std::pair<Point_3, Point_3>, EdgeWPtr> p2e;
      for (const EdgeSPtr& e : facet->edges()) {
        const Point_3& p0 = e->source()->point();
        const Point_3& p1 = e->target()->point();
        p2e.emplace((p0 < p1 ? std::make_pair(p0, p1) : std::make_pair(p1, p0)), e);
      }

      // Correspondence between the edges of the input mesh and the new edges in the polyhedron.
      // 'poly_vertices' is filled, starting at source() of the first edge, and
      // Facet::create() creates the i-th edge between vertices[i] and vertices[i+1]
      halfedge_descriptor h = halfedge(fd, tmesh), start_h = halfedge(fd, tmesh);
      do {
        const Point_3& p0 = get(CGAL::vertex_point, tmesh, source(h, tmesh));
        const Point_3& p1 = get(CGAL::vertex_point, tmesh, target(h, tmesh));
        EdgeWPtr ew = p2e.at((p0 < p1 ? std::make_pair(p0, p1) : std::make_pair(p1, p0)));
        e2e[edge(h, tmesh)] = ew;
        h = next(h, tmesh);
      } while (h != start_h);

      Vector_3 n = Polygon_mesh_processing::compute_face_normal(fd, tmesh);
      if (outward_offsetting)
        n = -n;

      Plane_3 pl (poly_vertices[0]->point(), n);
      facet->set_plane(pl);
      result->add_facet(facet);

      const FT weight = get(weight_pmap, fd);
      CGAL_assertion(weight > 0);

      SkelFacetDataSPtr data = Skeleton_facet_data::create(facet);
      data->set_speed(weight);
      data->set_input_face_id(facet_id_new);
    }

    for (const EdgeSPtr& edge : result->edges()) {
      if (!(edge->get_facet_L() && edge->get_facet_R())) {
        CGAL_SS3_IO_TRACE_V(1, "Warning: polyhedron has no closed boundary.");
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
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;
    CGAL::unordered_flat_map<edge_descriptor, EdgeWPtr> unused_e2e;
    CGAL::unordered_flat_map<face_descriptor, FacetWPtr> unused_f2f;
    return load(tmesh, unused_e2e, unused_f2f, np);
  }

  template <typename TriangleMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static PolyhedronSPtr convert(const TriangleMesh& tmesh,
                                const NamedParameters& np = CGAL::parameters::default_values())
  {
    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::get_parameter;

    CGAL_SS3_IO_TRACE_V(4, "Converting mesh...");

    using edge_descriptor = typename boost::graph_traits<TriangleMesh>::edge_descriptor;
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    const bool outward_offsetting = choose_parameter(get_parameter(np, internal_np::outward_offsetting), false);

    ConfigurationSPtr config = Configuration::get_instance();
    const bool merge_faces = config->get_Boolean("Preprocessing", "merge_coplanar_faces");
    const double epsilon = config->get_double("Preprocessing", "coplanarity_epsilon");

    CGAL::unordered_flat_map<edge_descriptor, EdgeWPtr> e2e;
    CGAL::unordered_flat_map<face_descriptor, FacetWPtr> f2f;

    PolyhedronSPtr polyhedron = load(tmesh, e2e, f2f, np);
    if (!merge_faces)
      return polyhedron;

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/coplanar_merge_before.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

#ifdef CGAL_SS3_DETECT_COPLANARITIES_WITH_NORMAL_CHANGE
    Transformation::merge_coplanar_facets(polyhedron);
#else // CGAL_SS3_DETECT_COPLANARITIES_WITH_NORMAL_CHANGE
    namespace PMP = CGAL::Polygon_mesh_processing;

    using halfedge_descriptor = typename boost::graph_traits<TriangleMesh>::halfedge_descriptor;
    using face_descriptor = typename boost::graph_traits<TriangleMesh>::face_descriptor;

    const CGAL::Bbox_3 bbox = PMP::bbox(tmesh);
    const FT diag_length = CGAL::approximate_sqrt(square(bbox.xmax() - bbox.xmin()) +
                                                  square(bbox.ymax() - bbox.ymin()) +
                                                  square(bbox.zmax() - bbox.zmin()));

    // Use shape detection to analyze the mesh
    std::vector<std::size_t> region_ids(num_faces(tmesh));
    boost::vector_property_map<Plane_3> plane_map; // supporting planes of the regions detected

    const FT max_distance = epsilon * diag_length;
    CGAL_SS3_IO_TRACE_V(8, "Region growing::max_distance = " << max_distance);

    // Detect planar regions in the mesh
    //
    // The cosine to '-1' is to ignore the angle change in consecutive faces and edges: elements
    // are part of the same region as long as they live within the same slab.
    // The reason behind this is that some input data (and results of this skeleton/offsetting
    // algorithm) can have nasty folds.

    PMP::region_growing_of_planes_on_faces(tmesh,
                                           CGAL::make_random_access_property_map(region_ids),
                                           CGAL::parameters::region_primitive_map(plane_map)
                                                            .maximum_distance(max_distance)
                                                            .cosine_of_maximum_angle(-1.));

    CGAL_SS3_IO_TRACE_CODE(for (face_descriptor f : faces(tmesh)))
    CGAL_SS3_IO_TRACE_V(16, "facet " << f << " is in region " << region_ids[f]);

#ifdef CGAL_SS3_DUMP_FILES
    static int region_dump_id = -1;
    utils::save_colored_mesh(tmesh, region_ids, "results/regions_" + std::to_string(++region_dump_id) + ".ply");
#endif // CGAL_SS3_DUMP_FILES

    // merge the facets incident to an unconstrained edge (i.e., the edge is interior to a region)
    std::vector<EdgeWPtr> edges_to_remove;
    for (edge_descriptor ed : edges(tmesh)) {
      halfedge_descriptor hd = halfedge(ed, tmesh);
      if (region_ids[face(hd, tmesh)] == region_ids[face(opposite(hd, tmesh), tmesh)]) {
        edges_to_remove.push_back(e2e.at(ed));
      }
    }

    Transformation::merge_facet_pairs(edges_to_remove, polyhedron);

    // merge facet pairs arbitrarily keeps one of the two facets, so we need the correct plane
    for (FacetSPtr f : polyhedron->facets()) {
      // abusing the fact that polyhedron and input mesh facets are in the same order...
      const Plane_3& pl = get(plane_map, region_ids[f->id()]);
      f->set_plane(outward_offsetting ? pl.opposite() : pl);
    }

    Transformation::sanitize(polyhedron); // remove degenerate vertices and facets
    CGAL_postcondition(polyhedron->is_consistent());

    polyhedron->initialize_all_IDs();

    // all vertices should be within the max distance of the chosen supporting plane
    for (const FacetSPtr& f : polyhedron->facets()) {
      for (const VertexSPtr& v : f->vertices()) {
        CGAL_postcondition(CGAL::squared_distance(f->get_plane(), v->point()) < CGAL::square(max_distance));
      }
    }

#endif // CGAL_SS3_DETECT_COPLANARITIES_WITH_NORMAL_CHANGE

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/coplanar_merge_after.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

#if 0
    // @todo test this again
    Transformation::truncate_precision(polyhedron);
#endif

#ifdef CGAL_SS3_DUMP_FILES
    IO::write_OBJ("results/converted.obj", polyhedron, parameters::do_not_triangulate_faces(true));
#endif

    return polyhedron;
  }

public:
  template <typename PolygonMesh,
            typename NamedParameters = CGAL::parameters::Default_named_parameters>
  static bool save(const PolyhedronSPtr& polyhedron,
                   PolygonMesh& pmesh,
                   const NamedParameters& np = CGAL::parameters::default_values())
  {
    namespace PMP = CGAL::Polygon_mesh_processing;

    using CGAL::parameters::choose_parameter;
    using CGAL::parameters::is_default_parameter;
    using CGAL::parameters::get_parameter;

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
    using Point = typename boost::property_traits<VPM>::value_type;

    CGAL_SS3_IO_TRACE_V(8, "Save polyhedron with " << polyhedron->vertices().size() << " vertices and "
                                                   << polyhedron->facets().size() << " facets");
    CGAL_precondition(polyhedron->is_consistent());

    // @todo do not systematically triangulate, but use this NP and if it is false,
    // only triangulate what is not representable otherwise (see code in PMP::remesh_planar_faces)
    // bool do_triangulate = !choose_parameter(get_parameter(np, CGAL::internal_np::do_not_triangulate_faces), false);

    std::vector<Point> points;
    std::vector<std::vector<std::size_t>> soup_faces;
    std::vector<FT> face_speeds;
    std::vector<std::size_t> face_input_ids;

    CGAL::unordered_flat_map<VertexSPtr, std::size_t> v_index_map;
    std::size_t vidx = 0;
    for (const VertexSPtr& vertex : polyhedron->vertices()) {
      points.push_back(vertex->point());
      v_index_map[vertex] = vidx++;
    }

    // Write facets
    auto weight_pmap = choose_parameter(get_parameter(np, CGAL::internal_np::face_weight),
                                        CGAL::Constant_property_map<std::size_t, FT>(1));

    std::map<face_descriptor, std::size_t> default_f2i;
    auto f2i = choose_parameter(get_parameter(np, internal_np::face_to_face_map),
                                boost::make_assoc_property_map(default_f2i));

    for (const FacetSPtr& facet : polyhedron->facets()) {
      CGAL_assertion(facet->id() != -1);
      CGAL_assertion(facet->edges().size() >= 3);

      FT speed = Hds_utils::get_speed(facet);
      std::size_t input_face_id = Hds_utils::get_input_face_id(facet);
      CGAL_SS3_IO_TRACE_V(32, "Saving facet " << facet->id() << " with speed " << speed << " and input face id " << input_face_id);

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
          CGAL_SS3_IO_TRACE_V(1, "Warning: degenerate edge at " << v0->point());
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
            CGAL_SS3_IO_TRACE_V(1, "Error: Intersection of constraints");
            CGAL_SS3_IO_TRACE_V(1, "While inserting " << v0->point() << " || " << v1->point());
            CGAL_SS3_IO_TRACE_V(1, facet->to_string());
            CGAL_assertion_msg(false, "Intersections in CDT2 are not allowed");
            return false;
          }
          ++ne;
        }
      }

      if(ne < 3) // degenerate facet
      {
        CGAL_SS3_IO_TRACE_V(1, "Warning: skipping degenerate facet");
        continue;
      }

      std::unordered_map<PCDT_FH, bool> in_domain_map;
      boost::associative_property_map<std::unordered_map<PCDT_FH, bool>> in_domain(in_domain_map);
      CGAL::mark_domain_in_triangulation(pcdt, in_domain);

      for (auto fh : pcdt.finite_face_handles()) {
        if(!get(in_domain, fh))
          continue;

        std::vector<std::size_t> tri(3);
        tri[0] = v_index_map[fh->vertex(0)->info()];
        tri[1] = v_index_map[fh->vertex(1)->info()];
        tri[2] = v_index_map[fh->vertex(2)->info()];
        soup_faces.emplace_back(std::move(tri));
        face_speeds.emplace_back(speed);
        face_input_ids.emplace_back(input_face_id);
      }
    }

    if(!PMP::is_polygon_soup_a_polygon_mesh(soup_faces))
    {
      CGAL_SS3_IO_TRACE("Warning: polygon soup does not describe a polygon mesh");
#ifdef CGAL_SS3_DUMP_FILES
      CGAL::IO::write_STL("results/nm_soup.stl", points, soup_faces);
#endif
      PMP::duplicate_non_manifold_edges_in_polygon_soup(points, soup_faces);
      CGAL_assertion(PMP::is_polygon_soup_a_polygon_mesh(soup_faces));
    }

    // Convert polygon soup to polygon mesh
    PMP::polygon_soup_to_polygon_mesh(points, soup_faces, pmesh, parameters::default_values(), np);

    // Transfer per-face properties in the same order as the soup faces were added
    std::size_t fi = 0;
    for (const face_descriptor& f : faces(pmesh)) {
      if (fi >= face_speeds.size()) break;
      put(weight_pmap, f, face_speeds[fi]);
      if constexpr (!is_default_parameter<NamedParameters, internal_np::face_to_face_map_t>::value)
        put(f2i, f, face_input_ids[fi]);
      ++fi;
    }

    return true;
  }
};

} // namespace IO
} // namespace Straight_skeletons_3
} // namespace CGAL

#endif /* CGAL_STRAIGHT_SKELETON_3_IO_SURFACE_MESH_IO_H */
