#include <CGAL/config.h>

// #define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#include <CGAL/assertions.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Conforming_Delaunay_triangulation_3.h>
#include <CGAL/Constrained_triangulation_3/internal/read_polygon_mesh_for_cdt_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/enum.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Random.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO/PLY.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/use.h>
#include <CGAL/utility.h>

#include <CGAL/bisect_failures.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>

#include <boost/graph/graph_traits.hpp>
#include <boost/iterator/function_output_iterator.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <ostream>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#if __has_include(<version>)
#  include <version>
#endif

#if CGAL_USE_CLI11
#  include <CLI/CLI.hpp>
#endif

#ifndef CGAL_CDT_3_USE_EPECK
#  define CGAL_CDT_3_USE_EPECK 0
#endif

namespace CGAL::CDT_3::test::cdt_3_from_off {

using K = std::conditional_t<CGAL_CDT_3_USE_EPECK == 0,
                             CGAL::Exact_predicates_inexact_constructions_kernel,
                             CGAL::Exact_predicates_exact_constructions_kernel>;

using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3<K>;
using Point = K::Point_3;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor   = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor   = boost::graph_traits<Mesh>::face_descriptor;

template <typename Rep, typename Period>
auto milliseconds(std::chrono::duration<Rep, Period> duration) {
  return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
};

template <typename Clock, typename Duration>
auto milliseconds_since(std::chrono::time_point<Clock, Duration> start_time) {
  return milliseconds(Clock::now() - start_time);
};

using clock = std::chrono::high_resolution_clock;

struct CDT_options
{
  int         verbose_level                       = 0;
  bool        quiet                               = false;
  bool        merge_facets                        = true;
  bool        merge_facets_old_method             = false;
  bool        reject_self_intersections           = false;
  bool        repair_mesh                         = true;
  bool        read_mesh_with_operator             = false;
  bool        debug_Steiner_points                = false;
  bool        debug_Steiner_points_construction   = false;
  unsigned    debug_move_Steiner_vertices         = 0;
  bool        debug_input_faces                   = false;
  bool        debug_missing_regions               = false;
  bool        debug_regions                       = false;
  bool        debug_copy_triangulation_into_hole  = false;
  bool        debug_validity                      = false;
  bool        debug_finite_edges_map              = false;
  bool        debug_subconstraints_to_conform     = false;
  bool        debug_verbose_special_cases         = false;
  bool        debug_encroaching_vertices          = false;
  bool        debug_conforming_validation         = false;
  bool        debug_constraint_hierarchy          = false;
  bool        debug_geometric_errors              = false;
  bool        debug_polygon_insertion             = false;
  bool        debug_restore_faces                 = false;
  bool        use_finite_edges_map                = false;
  bool        move_Steiner_vertices_to_the_volume = true;
  bool        allow_moving_Steiner_vertices_to_create_negative_tets = false;
  bool        use_epeck_for_normals               = false;
  bool        use_epeck_for_Steiner_points        = false;
  bool        call_is_valid                       = true;
  bool        bisect_failures                     = false;
  double      vertex_vertex_epsilon               = 0.; // 1e-14;
  double      segment_vertex_epsilon              = 0.; // 1e-14;
  double      coplanar_polygon_max_angle          = 5.1;
  double      coplanar_polygon_max_distance       = 1e-6;
  std::string input_filename                      = CGAL::data_file_path("meshes/mpi.off");
  std::string output_filename                     {"dump.off"};
  std::string dump_patches_after_merge_filename   {};
  std::string dump_surface_mesh_after_merge_filename{};
  std::string dump_patches_borders_prefix         {};
  std::string dump_after_conforming_filename      {};

  CDT_options(int argc, char* argv[]);
};

CDT_options::CDT_options(int argc, char* argv[]) {
#if CGAL_USE_CLI11
  CLI::App app{"Build and process a constrained Delaunay triangulation from an OFF mesh"};
  app.set_help_flag("-h,--help", "print this help");
  app.allow_extras(false);

  app.add_flag("-V,--verbose", verbose_level,
               "verbose (can be used several times)")
      ->multi_option_policy(CLI::MultiOptionPolicy::Sum);

  app.add_flag("--quiet", quiet, "do not print anything");

  app.add_flag("--merge-facets,!--no-merge-facets", merge_facets,
               "merge facets into patches");
  app.add_flag("--merge-facets-old", merge_facets_old_method,
               "merge facets using the old method");
  app.add_flag("--bisect", bisect_failures, "bisect failures");
  app.add_flag("--reject-self-intersections", reject_self_intersections,
               "reject self-intersecting polygon soups");
  app.add_flag("--read-mesh-with-operator", read_mesh_with_operator,
               "read the mesh with operator>>");
  app.add_flag("--repair,!--no-repair", repair_mesh,
               "repair the mesh while reading polygon soup");
  app.add_flag("--call-is-valid,!--no-is-valid", call_is_valid,
               "call is_valid checks");

  app.add_option("--vertex-vertex-epsilon", vertex_vertex_epsilon,
                 "epsilon for vertex-vertex min distance");
  app.add_option("--segment-vertex-epsilon", segment_vertex_epsilon,
                 "epsilon for segment-vertex min distance");
  app.add_option("--coplanar-polygon-max-angle", coplanar_polygon_max_angle,
                 "max angle for coplanar polygons");
  app.add_option("--coplanar-polygon-max-distance", coplanar_polygon_max_distance,
                 "max distance for coplanar polygons");

  app.add_option("--dump-patches-after-merge", dump_patches_after_merge_filename,
                 "dump patches after merging facets in PLY");
  app.add_option("--dump-surface-mesh-after-merge", dump_surface_mesh_after_merge_filename,
                 "dump surface mesh after merging facets in OFF");
  app.add_option("--dump-patches-borders-prefix", dump_patches_borders_prefix,
                 "dump patches borders");
  app.add_option("--dump-after-conforming", dump_after_conforming_filename,
                 "dump mesh after conforming in OFF");

  app.add_flag("--debug-Steiner-points", debug_Steiner_points,
               "debug Steiner point insertion");
  app.add_flag("--debug-Steiner-points-construction", debug_Steiner_points_construction,
               "debug Steiner point construction");
  auto debug_move_Steiner_vertices_option =
      app.add_option("--debug-move-Steiner-vertices-to-the-volume", debug_move_Steiner_vertices,
                     "debug moving Steiner vertices to the volume (level 0..3)")
          ->check(CLI::Range(0, 3))
          ->expected(0, 1);
  app.add_flag("--debug-input-faces", debug_input_faces, "debug input faces");
  app.add_flag("--debug-missing-regions", debug_missing_regions, "debug missing regions");
  app.add_flag("--debug-regions", debug_regions, "debug regions");
  app.add_flag("--debug_copy_triangulation_into_hole", debug_copy_triangulation_into_hole,
               "debug copy_triangulation_into_hole");
  app.add_flag("--debug-validity", debug_validity,
               "add is_valid checks after modifications to the TDS");
  app.add_flag("--debug-finite-edges-map", debug_finite_edges_map,
               "debug the use of a hash map for finite edges");
  app.add_flag("--debug-subconstraints-to-conform", debug_subconstraints_to_conform,
               "debug subconstraints to conform");
  app.add_flag("--debug-verbose-special-cases", debug_verbose_special_cases,
               "debug verbose output for special cases");
  app.add_flag("--debug-encroaching-vertices", debug_encroaching_vertices,
               "debug encroaching vertices computation");
  app.add_flag("--debug-conforming-validation", debug_conforming_validation,
               "debug edge conforming validation");
  app.add_flag("--debug-constraint-hierarchy", debug_constraint_hierarchy,
               "debug constraint hierarchy operations");
  app.add_flag("--debug-geometric-errors", debug_geometric_errors,
               "debug geometric error handling");
  app.add_flag("--debug-polygon-insertion", debug_polygon_insertion,
               "debug polygon insertion process");
  app.add_flag("--debug-restore-faces", debug_restore_faces,
               "debug face restoration process");

  app.add_flag("--use-finite-edges-map", use_finite_edges_map,
               "use a hash map for finite edges");
  app.add_flag("--move-Steiner-vertices-to-the-volume,!--no-move-Steiner-vertices-to-the-volume",
               move_Steiner_vertices_to_the_volume,
               "move Steiner vertices to the volume");
  app.add_flag("--allow-moving-Steiner-vertices-to-create-negative-tets",
               allow_moving_Steiner_vertices_to_create_negative_tets,
               "allow moving Steiner vertices to create negative volume tetrahedra");
  app.add_flag("--use-epeck-for-normals,!--no-use-epeck-for-normals", use_epeck_for_normals,
               "use exact kernel for normal computations");
  app.add_flag("--use-epeck-for-Steiner-points,!--no-use-epeck-for-Steiner-points",
               use_epeck_for_Steiner_points,
               "use exact kernel for Steiner point computations");

  app.add_option("input.off", input_filename, "input mesh");
  app.add_option("output.off", output_filename, "output mesh");

  for(auto* option : app.get_options()) {
    option->always_capture_default(true);
  }

  app.get_formatter()->enable_default_flag_values(true);

  try {
    app.parse(argc, argv);
  } catch(const CLI::ParseError& e) {
    std::exit(app.exit(e));
  }

  if(*debug_move_Steiner_vertices_option && debug_move_Steiner_vertices == 0) {
    debug_move_Steiner_vertices = 1;
  }

  if(merge_facets_old_method) {
    merge_facets = true;
  }
#endif // CGAL_USE_CLI11
}

CGAL::CDT_3::Debug_options cdt_debug_options(const CDT_options& options) {
  CGAL::CDT_3::Debug_options cdt_debug;
  cdt_debug.Steiner_points(options.debug_Steiner_points || options.verbose_level > 0);
  cdt_debug.Steiner_points_construction(options.debug_Steiner_points_construction);
  cdt_debug.input_faces(options.debug_input_faces);
  cdt_debug.missing_region(options.verbose_level > 1 || options.debug_missing_regions);
  cdt_debug.regions(options.debug_regions);
  cdt_debug.validity(options.debug_validity);
  cdt_debug.finite_edges_map(options.debug_finite_edges_map);
  cdt_debug.subconstraints_to_conform(options.debug_subconstraints_to_conform);
  cdt_debug.verbose_special_cases(options.debug_verbose_special_cases);
  cdt_debug.encroaching_vertices(options.debug_encroaching_vertices);
  cdt_debug.conforming_validation(options.debug_conforming_validation);
  cdt_debug.move_Steiner_vertices(options.debug_move_Steiner_vertices);
  cdt_debug.move_Steiner_vertices_allow_negative_tets(options.allow_moving_Steiner_vertices_to_create_negative_tets);
  cdt_debug.constraint_hierarchy(options.debug_constraint_hierarchy);
  cdt_debug.geometric_errors(options.debug_geometric_errors);
  cdt_debug.polygon_insertion(options.debug_polygon_insertion);
  cdt_debug.restore_faces(options.debug_restore_faces);
  cdt_debug.copy_triangulation_into_hole(options.debug_copy_triangulation_into_hole);
  cdt_debug.use_finite_edges_map(options.use_finite_edges_map);
  cdt_debug.display_statistics(!options.quiet);
  cdt_debug.use_epeck_for_normals(options.use_epeck_for_normals);
  cdt_debug.use_epeck_for_Steiner_points(options.use_epeck_for_Steiner_points);
  cdt_debug.set_segment_vertex_epsilon(options.segment_vertex_epsilon);
  cdt_debug.set_vertex_vertex_epsilon(options.vertex_vertex_epsilon);

  return cdt_debug;
}

auto compute_bounding_box(const Mesh& mesh, const CDT_options& options) {
  const auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  double d_x = bbox.xmax() - bbox.xmin();
  double d_y = bbox.ymax() - bbox.ymin();
  double d_z = bbox.zmax() - bbox.zmin();

  const double bbox_max_width = (std::max)(d_x, (std::max)(d_y, d_z));

  if(!options.quiet) {
    double epsilon = options.vertex_vertex_epsilon;
    std::cout << "Bbox width                   : " << bbox_max_width << '\n'
              << "Epsilon                      : " << epsilon << '\n'
              << "Epsilon * Bbox width         : " << epsilon * bbox_max_width << "\n\n";
  }

  return bbox;
}

std::function<void()> create_output_finalizer(const CDT& cdt, const CDT_options& options) {
  return [&cdt, &options]() {
    auto _ = CGAL::CDT_3_OUTPUT_TASK_guard();
    {
      if(std::uncaught_exceptions() > 0) {
        std::cerr << "create_output_finalizer: an exception is being thrown, the triangulation may be left in an "
                     "invalid state.\n";
        return;
      }
      auto dump_to_medit = [](std::string fname,
                              const std::vector<K::Point_3> &points,
                              const std::vector<std::array<std::size_t, 4>> &indexed_triangles,
                              const std::vector<std::array<std::size_t, 4>> &indexed_tetra,
                              const std::vector<std::size_t> &cell_ids)
      {
        std::ofstream out(fname);
        out.precision(17);
        out << "MeshVersionFormatted 1\nDimension 3\nVertices\n";
        out << points.size() << "\n";
        for (const K::Point_3& p : points)
          out << p << " 0\n";
        out << "Triangles\n" << indexed_triangles.size() << "\n";
        for (const auto& tri : indexed_triangles)
          out << tri[0]+1 << " " << tri[1]+1 << " " << tri[2]+1 << " " << tri[3] << "\n";
        out << "Tetrahedra\n";
        out << indexed_tetra.size() << "\n";
        for (std::size_t k=0;k<indexed_tetra.size(); ++k)
          out << indexed_tetra[k][0]+1 << " "
              << indexed_tetra[k][1]+1 << " "
              << indexed_tetra[k][2]+1 << " "
              << indexed_tetra[k][3]+1 << " " << cell_ids[k] << "\n";
        out <<"End\n";
      };

      const auto& tr = cdt.triangulation();
      using Tr = CDT::Triangulation;

      std::unordered_map<Tr::Cell_handle, int /*Subdomain_index*/> cells_map;
      for(auto ch : tr.all_cell_handles())
      {
        cells_map[ch] = 1;
      }

      std::stack<decltype(tr.infinite_cell())> stack;
      stack.push(tr.infinite_cell());
      while (!stack.empty())
      {
        auto ch = stack.top();
        stack.pop();
        cells_map[ch] = 0;
        for (int i = 0; i < 4; ++i)
        {
          if(ch->ccdt_3_data().is_facet_constrained(i))
            continue;
          auto n = ch->neighbor(i);
          if (cells_map[n] == 1)
          {
            stack.push(n);
          }
        }
      }

      const bool export_all_the_bbox =
          std::all_of(cells_map.begin(), cells_map.end(), [](const auto& kv) { return kv.second == 0; });
      if (export_all_the_bbox) {
        for(auto ch : tr.finite_cell_handles())
        {
          cells_map[ch] = 1;
        }
      }

      std::vector<K::Point_3> points(cdt.number_of_vertices());
      std::size_t idx = 0;
      for(auto v: cdt.finite_vertex_handles()) {
        // renumber the vertices: there might be holes in the numbering due to removed vertices
        v->set_time_stamp(++idx);
        points.at(idx - 1) = v->point();
      }
      std::vector<std::array<std::size_t, 4>> indexed_triangles;
      std::size_t max_polygon_id = 0;
      for(auto f : cdt.constrained_facets()) {
        const auto verts = tr.vertices(f);
        const auto polygon_id = static_cast<std::size_t>(cdt.face_constraint_index(f));
        indexed_triangles.push_back({verts[0]->time_stamp() - 1,
                                     verts[1]->time_stamp() - 1,
                                     verts[2]->time_stamp() - 1,
                                     polygon_id});
        if (polygon_id > max_polygon_id) {
          max_polygon_id = polygon_id;
        }
      }
      if(export_all_the_bbox) {
        auto inf_v = tr.infinite_vertex();
        tr.incident_cells(inf_v, boost::make_function_output_iterator([&](Tr::Cell_handle ch) {
                            auto facet_index = ch->index(inf_v);
                            const auto verts = tr.vertices(Tr::Facet{ch, facet_index});
                            indexed_triangles.push_back({verts[0]->time_stamp() - 1,
                                                         verts[1]->time_stamp() - 1,
                                                         verts[2]->time_stamp() - 1,
                                                         max_polygon_id + 1 /*a new polygon id for the bbox*/});
                          }));
      }
      std::vector<std::array<std::size_t, 4>> indexed_tetra;
      indexed_tetra.reserve(cdt.number_of_cells());
      for(auto ch: cdt.finite_cell_handles()) {
        if(cells_map[ch] > 0) {
          indexed_tetra.push_back({ch->vertex(0)->time_stamp() - 1,
                                   ch->vertex(1)->time_stamp() - 1,
                                   ch->vertex(2)->time_stamp() - 1,
                                   ch->vertex(3)->time_stamp() - 1});
        }
      }
      std::vector<std::size_t> cell_ids(indexed_tetra.size(), 1);
      dump_to_medit(options.output_filename + ".mesh", points, indexed_triangles, indexed_tetra, cell_ids);
      cdt.dump_constrained_facets_to_off(options.output_filename);
    }
    {
      cdt.save_binary_file(options.output_filename + ".binary.cgal");
    }
  };
}


template<typename PatchIdMap, typename VertexSelectedMap, typename EdgeBorderMap, typename VertexPointMap>
struct Mesh_property_maps {
  PatchIdMap patch_id_map;
  VertexSelectedMap v_selected_map;
  EdgeBorderMap edge_is_border_of_patch_map;
  VertexPointMap mesh_vertex_point_map;
};

// CTAD deduction guide
template<typename PatchIdMap, typename VertexSelectedMap, typename EdgeBorderMap, typename VertexPointMap>
Mesh_property_maps(PatchIdMap, VertexSelectedMap, EdgeBorderMap, VertexPointMap)
        -> Mesh_property_maps<PatchIdMap, VertexSelectedMap, EdgeBorderMap, VertexPointMap>;

auto setup_mesh_property_maps(Mesh& mesh) {
  auto [patch_id_map, patch_id_map_ok] = mesh.add_property_map<face_descriptor, int>("f:patch_id", -2);
  assert(patch_id_map_ok); CGAL_USE(patch_id_map_ok);
  auto [v_selected_map, v_selected_map_ok] = mesh.add_property_map<vertex_descriptor, bool>("v:selected", false);
  assert(v_selected_map_ok); CGAL_USE(v_selected_map_ok);
  auto [edge_is_border_of_patch_map, edge_is_border_of_patch_map_ok] =
      mesh.add_property_map<edge_descriptor, bool>("e:is_border_of_patch", false);
  assert(edge_is_border_of_patch_map_ok); CGAL_USE(edge_is_border_of_patch_map_ok);
  auto mesh_vertex_point_map = get(CGAL::vertex_point, mesh);

  return Mesh_property_maps{patch_id_map, v_selected_map, edge_is_border_of_patch_map, mesh_vertex_point_map};
}

using Borders_of_patches = std::vector<std::vector<std::pair<vertex_descriptor, vertex_descriptor>>>;

template<typename MeshPropertyMaps>
auto extract_patch_edges(Mesh& mesh, MeshPropertyMaps pmaps, int number_of_patches) {
  Borders_of_patches patch_edges(number_of_patches);

  for(auto h: halfedges(mesh))
  {
    if(is_border(h, mesh)) continue;
    auto f = face(h, mesh);
    auto patch_id = get(pmaps.patch_id_map, f);
    auto opp = opposite(h, mesh);
    if(is_border(opp, mesh) || patch_id != get(pmaps.patch_id_map, face(opp, mesh))) {
      auto va = source(h, mesh);
      auto vb = target(h, mesh);
      patch_edges[patch_id].emplace_back(va, vb);
      put(pmaps.v_selected_map, va, true);
      put(pmaps.v_selected_map, vb, true);
    }
  }

  return patch_edges;
}

template <typename MeshPropertyMaps>
int merge_facets_region_growing(Mesh& mesh,
                                MeshPropertyMaps pmaps,
                                double coplanar_polygon_max_distance,
                                double coplanar_polygon_max_angle,
                                const std::string& dump_surface_mesh_after_merge_filename) {
  namespace np = CGAL::parameters;
  int number_of_patches = static_cast<int>(CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces(
      mesh, pmaps.patch_id_map,
      np::maximum_distance(coplanar_polygon_max_distance)
         .maximum_angle(coplanar_polygon_max_angle)));
  for(auto f: faces(mesh)) {
    if(get(pmaps.patch_id_map, f) < 0) {
      std::cerr << "warning: face " << f << " has no patch id! Reassign it to " << number_of_patches << '\n';
      for(auto h: CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
        std::cerr << "  " << target(h, mesh) << ", point " << mesh.point(target(h, mesh)) << '\n';
      }
      put(pmaps.patch_id_map, f, number_of_patches++);
    }
  }
  if(!dump_surface_mesh_after_merge_filename.empty()) {
    static constexpr auto max_size_t = (std::numeric_limits<std::size_t>::max)();
    auto [corner_id_map, corner_id_map_ok] =
        mesh.add_property_map<vertex_descriptor, std::size_t>("v:corner_id", max_size_t);
    assert(corner_id_map_ok);
    CGAL_USE(corner_id_map_ok);
    const auto nb_corners = CGAL::Polygon_mesh_processing::detect_corners_of_regions(
        mesh, pmaps.patch_id_map, number_of_patches, corner_id_map,
        np::maximum_distance(coplanar_polygon_max_distance)
            .maximum_angle(coplanar_polygon_max_angle)
            .edge_is_constrained_map(pmaps.edge_is_border_of_patch_map));
    Mesh merged_mesh;
    CGAL::Polygon_mesh_processing::remesh_almost_planar_patches(
        mesh, merged_mesh, number_of_patches, nb_corners, pmaps.patch_id_map,
        corner_id_map, pmaps.edge_is_border_of_patch_map,
        CGAL::parameters::default_values(),
        CGAL::parameters::do_not_triangulate_faces(true));
    mesh.remove_property_map(corner_id_map);
    std::ofstream out(dump_surface_mesh_after_merge_filename);
    out.precision(17);
    out << merged_mesh;
  }

  return number_of_patches;
}

template<typename MeshPropertyMaps>
int merge_facets_old_method(Mesh& mesh, MeshPropertyMaps pmaps, int initial_number_of_patches) {
  int number_of_patches = initial_number_of_patches;

  for(auto f: faces(mesh))
  {
    if(get(pmaps.patch_id_map, f) >= 0) continue;
    std::stack<face_descriptor> f_stack;
    f_stack.push(f);
    while(!f_stack.empty()) {
      auto f = f_stack.top();
      f_stack.pop();
      if(get(pmaps.patch_id_map, f) >= 0) continue;
      put(pmaps.patch_id_map, f, number_of_patches);
      for(auto h: CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
        auto opp = opposite(h, mesh);
        if(is_border_edge(opp, mesh)) {
          continue;
        }
        auto n = face(opp, mesh);
        auto a = get(pmaps.mesh_vertex_point_map, source(h, mesh));
        auto b = get(pmaps.mesh_vertex_point_map, target(h, mesh));
        auto c = get(pmaps.mesh_vertex_point_map, target(next(h, mesh), mesh));
        auto d = get(pmaps.mesh_vertex_point_map, target(next(opp, mesh), mesh));
        if(CGAL::orientation(a, b, c, d) != CGAL::COPLANAR) {
          continue;
        }
        if(get(pmaps.patch_id_map, n) >= 0) continue;
        f_stack.push(n);
      }
    }
    ++number_of_patches;
  }

  return number_of_patches;
}

template<typename MeshPropertyMaps>
Borders_of_patches maybe_merge_facets(
    Mesh& mesh,
    const CDT_options& options,
    MeshPropertyMaps pmaps,
    double bbox_max_span)
{
  int number_of_patches = 0;
  Borders_of_patches patch_edges;

  if(!options.merge_facets) return patch_edges;

  {
    auto _ = CGAL::CDT_3_MERGE_FACETS_TASK_guard();
    auto start_time = clock::now();

    if(options.merge_facets_old_method) {
      number_of_patches = merge_facets_old_method(mesh, pmaps, number_of_patches);
    } else {
      number_of_patches = merge_facets_region_growing(
          mesh, pmaps, options.coplanar_polygon_max_distance * bbox_max_span,
          options.coplanar_polygon_max_angle, options.dump_surface_mesh_after_merge_filename);
    }
    patch_edges = extract_patch_edges(mesh, pmaps, number_of_patches);
    if (!options.quiet) {
      std::cout << "[timings] found " << number_of_patches << " patches in "
                << milliseconds_since(start_time) << " ms\n";
    }
  }

  if(options.dump_patches_after_merge_filename.empty()) return patch_edges;

  {
    auto _ = CGAL::CDT_3_OUTPUT_TASK_guard();
    std::ofstream out(options.dump_patches_after_merge_filename);
    CGAL::IO::write_PLY(out, mesh, CGAL::parameters::stream_precision(17));
  }

  return patch_edges;
}

template<typename BordersOfPatches, typename MeshPropertyMaps>
void dump_patches_borders(const BordersOfPatches& patch_edges,
                          const MeshPropertyMaps& pmaps,
                          const std::string& dump_patches_borders_prefix) {
  std::set<std::pair<vertex_descriptor, vertex_descriptor>> all_edges;
  for(auto i = 0u; i < patch_edges.size(); ++i) {
    std::stringstream ss;
    ss << dump_patches_borders_prefix << i << ".polylines.txt";
    std::ofstream out(ss.str());
    out.precision(17);
    const auto& edges = patch_edges[i];
    for(auto [va, vb]: edges) {
      all_edges.insert(CGAL::make_sorted_pair(va, vb));
    }
    std::cerr << "Patch p#" << i << " has " << edges.size() << " edges\n";
    const auto polylines = CGAL::segment_soup_to_polylines(edges);
    for(const auto& polyline: polylines) {
      out << polyline.size() << "    ";
      for(auto v: polyline) {
        out << get(pmaps.mesh_vertex_point_map, v) << "  ";
      }
      out << '\n';
    }
    out.close();
    std::cerr << "  " << polylines.size() << " polylines\n";
    for(const auto& polyline: polylines) {
      std::cerr << "    - " << polyline.size() << " vertices\n";
      assert(polyline.front() == polyline.back());
    }
  }
  std::stringstream ss;
  ss << dump_patches_borders_prefix << "all_edges.polylines.txt";
  std::ofstream out(ss.str());
  out.precision(17);
  const auto polylines = CGAL::segment_soup_to_polylines(all_edges);
  for(const auto& polyline: polylines) {
    out << polyline.size() << "    ";
    for(auto v: polyline) {
      out << get(pmaps.mesh_vertex_point_map, v) << "  ";
    }
    out << '\n';
  }
}

int go(Mesh mesh, CDT_options options) {
  CGAL::CDT_3::Debug_options cdt_debug = cdt_debug_options(options);

  // Compute bbox (used for distance checks scaling)
  auto bbox = compute_bounding_box(mesh, options);
  auto bbox_max_span = (std::max)(bbox.x_span(), (std::max)(bbox.y_span(), bbox.z_span()));

  auto pmaps = setup_mesh_property_maps(mesh);

  auto patch_edges = maybe_merge_facets(mesh, options, pmaps, bbox_max_span);

  if(!options.dump_patches_borders_prefix.empty()) {
    auto _ = CGAL::CDT_3_OUTPUT_TASK_guard();
    dump_patches_borders(patch_edges, pmaps, options.dump_patches_borders_prefix);
  }

  auto dump_mesh_with_steiner_points = [&](const auto& cdt) {
    if(options.dump_after_conforming_filename.empty()) {
      return;
    }
    auto _ = CGAL::CDT_3_OUTPUT_TASK_guard();
    using Vertex_index = Mesh::Vertex_index;
    [[maybe_unused]] std::size_t time_stamp_counter = 0u;
    for(auto v : cdt.finite_vertex_handles()) {
      [[maybe_unused]] const auto time_stamp = v->time_stamp();
      assert(++time_stamp_counter == time_stamp);
      if(!v->ccdt_3_data().is_Steiner_vertex_on_edge())
        continue;
      const auto [va, vb] = cdt.ancestors_of_Steiner_vertex_on_edge(v);
      const auto index_va = Vertex_index{static_cast<unsigned>(va->time_stamp() - 1)};
      const auto index_vb = Vertex_index{static_cast<unsigned>(vb->time_stamp() - 1)};
      auto [it, end] = CGAL::halfedges_around_source(index_va, mesh);
      it = std::find_if(it, end, [&mesh, index_vb](auto he) { return target(he, mesh) == index_vb; });
      CGAL_assertion(it != end);
      auto he = CGAL::Euler::split_edge(*it, mesh);
      auto mesh_v = target(he, mesh);
      put(pmaps.mesh_vertex_point_map, mesh_v, v->point());
      assert(mesh_v == Vertex_index{static_cast<unsigned>(time_stamp - 1)});
    }
    std::ofstream out_mesh(options.dump_after_conforming_filename);
    out_mesh.precision(17);
    out_mesh << mesh;
    out_mesh.close();
  };

  // Build CDT directly from the polygon mesh using the official API
  CDT cdt;
  auto output_on_exit_scope_guard = CGAL::make_scope_exit(create_output_finalizer(cdt, options));
  auto build_start = clock::now();
  cdt = CDT{std::invoke([&]() {
    auto np = CGAL::parameters::debug(cdt_debug).visitor(std::ref(dump_mesh_with_steiner_points));
    return options.merge_facets ? CDT(mesh, np.plc_face_id(pmaps.patch_id_map)) : CDT(mesh, np);
  })};
  if(!options.quiet) {
    std::cout << "[timings] built CDT from mesh in "
              << milliseconds_since(build_start) << " ms\n";
    std::cout << cdt.statistics() << "\n\n";
  }

  // Validation: check conforming status and validity
  {
    auto _ = CGAL::CDT_3_VALIDATION_TASK_guard();
    CGAL_assertion(!options.call_is_valid || cdt.is_valid(true));
    CGAL_assertion(!options.call_is_valid || cdt.is_conforming());
  }

  // Move Steiner vertices to the volume if requested
  if(options.move_Steiner_vertices_to_the_volume) {
    auto move_vertices_guard = CGAL::CDT_3_MOVE_STEINER_VERTICES_TASK_guard();
    cdt.move_Steiner_vertices_to_the_volume();
    if(!options.quiet) {
      std::cout << "\n[timings] moved Steiner vertices to the volume in "
                    << move_vertices_guard.time_ms() << " ms\n";
      std::cout << cdt.statistics() << "\n";
    }
  }

  return EXIT_SUCCESS;
}

int bisect_errors(Mesh mesh, CDT_options options) {
  // Lambda to get the size (number of faces) of the mesh
  auto get_size = [](const Mesh& m) -> std::size_t {
    return m.number_of_faces();
  };

  // Lambda to simplify the mesh by removing faces in the range [start, end)
  auto simplify = [](Mesh& m, std::size_t start, std::size_t end) -> bool {
    // Remove faces from end-1 down to start (reverse order to maintain indices)
    auto face_begin = m.faces().begin();
    std::advance(face_begin, start);
    auto face_end = face_begin;
    std::advance(face_end, end - start);

    for(auto it = face_end; it != face_begin; ) {
      --it;
      CGAL::Euler::remove_face(halfedge(*it, m), m);
    }

    std::stringstream ss;
    ss.precision(17);
    ss << m;
    ss.seekg(0);
    Mesh simplified_mesh;
    ss >> simplified_mesh;
    m = std::move(simplified_mesh);

    return m.is_valid(true);
    return true;
  };

  // Lambda to run the test on the mesh
  auto run_test = [&options](const Mesh& m) -> int {
    return go(m, options);
  };

  // Lambda to save the mesh to a file
  auto save_mesh = [](const Mesh& m, const std::string& prefix) {
    static int count = 0;
    auto filename = prefix + "_mesh";
    if(prefix == "bad") {
      ++count;
      filename += "_";
      filename += std::to_string(count) + ".off";
    } else {
      filename += ".off";
    }
    std::ofstream out(filename);
    out.precision(17);
    out << m;
  };

  return CGAL::bisect_failures(mesh, get_size, simplify, run_test, save_mesh);
}

} // end of namespaces CGAL::CDT_3::test::cdt_3_from_off

int main(int argc, char* argv[]) {
  using namespace CGAL::CDT_3::test::cdt_3_from_off;
  std::cerr.precision(17);
  std::cout.precision(17);
  CGAL::get_default_random() = CGAL::Random(42);

  CDT_options options(argc, argv);

  CGAL::CDT_3_read_polygon_mesh_output<Mesh> read_mesh_result;
  auto start_time = clock::now();
  {
    auto _ = CGAL::CDT_3_READ_INPUT_TASK_guard();

    if(options.read_mesh_with_operator) {
      std::ifstream in(options.input_filename);
      if(!in) {
        std::cerr << "Cannot open file: " << options.input_filename << std::endl;
        return EXIT_FAILURE;
      }
      in >> *read_mesh_result.polygon_mesh;
      if(!in) {
        std::cerr << "Error reading mesh with operator>>" << std::endl;
        return EXIT_FAILURE;
      }
    } else {
      auto read_options = CGAL::parameters::repair_polygon_soup(options.repair_mesh).verbose(options.verbose_level);
      read_mesh_result = CGAL::read_polygon_mesh_for_cdt_3<Mesh>(options.input_filename, read_options);
    }

    if (!read_mesh_result.polygon_mesh)
    {
      std::cerr << "Not a valid input file." << std::endl;
      std::cerr << "Details:\n" << read_mesh_result.polygon_mesh.error() << std::endl;
      return EXIT_FAILURE;
    }
    if(!options.quiet) {
      std::cout << "[timings] read mesh \"" << options.input_filename << "\" in "
                << milliseconds_since(start_time)
                << " ms\n";
      std::cout << "Number of vertices: " << read_mesh_result.polygon_mesh->number_of_vertices() << '\n';
      std::cout << "Number of edges: " << read_mesh_result.polygon_mesh->number_of_edges() << '\n';
      std::cout << "Number of faces: " << read_mesh_result.polygon_mesh->number_of_faces() << "\n\n";

      if(!options.read_mesh_with_operator) {
        std::cout << "Processing was successful.\n";
        std::cout << "  Number of duplicated points: " << read_mesh_result.nb_of_duplicated_points << '\n';
        std::cout << "  Number of simplified polygons: " << read_mesh_result.nb_of_simplified_polygons << '\n';
        std::cout << "  Number of new polygons: " << read_mesh_result.nb_of_new_polygons << '\n';
        std::cout << "  Number of removed invalid polygons: " << read_mesh_result.nb_of_removed_invalid_polygons << '\n';
        std::cout << "  Number of removed duplicated polygons: " << read_mesh_result.nb_of_removed_duplicated_polygons << '\n';
        std::cout << "  Number of removed isolated points: " << read_mesh_result.nb_of_removed_isolated_points << '\n';
        std::cout << "  Polygon soup self-intersects: " << (read_mesh_result.polygon_soup_self_intersects ? "YES" : "no") << '\n';
        std::cout << "  Polygon mesh is manifold: " << (read_mesh_result.polygon_mesh_is_manifold ? "yes" : "NO") << '\n';
        std::cout << std::endl;
      }
    }
  }

  if(options.reject_self_intersections && read_mesh_result.polygon_soup_self_intersects) {
    std::cerr << "ERROR: input mesh self-intersects\n";
    return EXIT_FAILURE;
  }

  if(options.bisect_failures) {
    return bisect_errors(std::move(*read_mesh_result.polygon_mesh), options);
  }

  auto exit_code = go(std::move(*read_mesh_result.polygon_mesh), options);
  if(!options.quiet) {
    std::cout << "[timings] total time: " << milliseconds_since(start_time) << " ms\n";
    if(exit_code != 0) std::cout << "ERROR with exit code " << exit_code << '\n';
  }
  return exit_code;
}
