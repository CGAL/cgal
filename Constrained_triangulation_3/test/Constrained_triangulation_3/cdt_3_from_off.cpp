#ifdef _MSC_VER
#pragma warning(disable: 4455)
#endif

#if defined(CGAL_DEBUG_CDT_3) && !__has_include(<format>)
#undef CGAL_DEBUG_CDT_3
#endif

// #define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#define NO_TRY_CATCH 1
// #define CGAL_DEBUG_CDT_3 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_vertex_base_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Constrained_triangulation_3/internal/read_polygon_mesh_for_cdt_3.h>
#include <CGAL/IO/File_binary_mesh_3.h>
#include <CGAL/utility.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>

#include <boost/graph/graph_traits.hpp>

#include <cassert>
#include <chrono>
#include <exception>
#include <fstream>
#include <optional>
#include <string_view>
#include <string>
#include <utility>
#include <vector>

#if CGAL_CXX20 && __cpp_lib_concepts >= 201806L && __cpp_lib_ranges >= 201911L
#include <ranges>
#endif

#if CGAL_CDT_3_USE_EPECK

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;

#else // use Epick

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

#endif // use Epick

struct Vb : public CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K> {};
struct Cb : public CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K> {};
struct Tds: public CGAL::Triangulation_data_structure_3<Vb, Cb> {};
using Base_triantulation = CGAL::Delaunay_triangulation_3<K, Tds>;
using CDT = CGAL::Conforming_constrained_Delaunay_triangulation_3_impl<Base_triantulation>;
using Point = Base_triantulation::Point;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor   = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor   = boost::graph_traits<Mesh>::face_descriptor;

void help(std::ostream& out) {
  out << R"!!!!!(
Usage: cdt_3_from_off [options] input.off output.off

  input.off: input mesh
  output.off: output mesh

  --merge-facets/--no-merge-facets: merge facets into patches (set by default)
  --merge-facets-old: merge facets using the old method
  --use-new-cavity-algorithm/--use-old-cavity-algorithm: use new or old cavity algorithm (default: new)
  --failure-expression <expression>: expression to detect bad meshratio)
  --ratio <double>: ratio of faces to remove (default: 0.1), if --failure-expression is used
  --vertex-vertex-epsilon <double>: epsilon for vertex-vertex min distance (default: 1e-6)
  --segment-vertex-epsilon <double>: epsilon for segment-vertex min distance (default: 0)
  --coplanar-polygon-max-angle <double>: max angle for coplanar polygons (default: 1)
  --coplanar-polygon-max-distance <double>: max distance for coplanar polygons (default: 1e-6)

  --dump-patches-after-merge <filename.ply>: dump patches after merging facets in PLY
  --dump-surface-mesh-after-merge <filename.off>: dump surface mesh after merging facets in OFF
  --dump-patches-borders-prefix <filenames_prefix>: dump patches borders
  --dump-after-conforming <filename.off>: dump mesh after conforming in OFF

  --no-repair: do not repair the mesh
  --read-mesh-with-operator: read the mesh with operator>>
  --reject-self-intersections: reject self-intersecting polygon soups
  --no-is-valid: do not call is_valid checks
  --debug-input-faces: debug input faces
  --debug-missing-regions: debug missing regions
  --debug-regions: debug regions
  --debug_copy_triangulation_into_hole: debug copy_triangulation_into_hole
  --debug-validity: add is_valid checks after modifications to the TDS
  --debug-finite-edges-map: debug the use of a hash map for finite edges
  --debug-subconstraints-to-conform: debug subconstraints to conform
  --debug-verbose-special-cases: debug verbose output for special cases
  --debug-encroaching-vertices: debug encroaching vertices computation
  --debug-conforming-validation: debug edge conforming validation
  --debug-constraint-hierarchy: debug constraint hierarchy operations
  --debug-geometric-errors: debug geometric error handling
  --debug-polygon-insertion: debug polygon insertion process
  --use-finite-edges-map: use a hash map for finite edges (default: false)

  --verbose/-V: verbose (can be used several times)
  --quiet: do not print anything
  --help/-h: print this help
)!!!!!";
}

[[noreturn]] void error(std::string_view message, std::string_view extra = "") {
  std::cerr << "Error: " << message << extra << '\n';
  help(std::cerr);
  std::exit(EXIT_FAILURE);
}

struct CDT_options
{
  int         verbose_level                       = 0;
  bool        need_help                           = false;
  bool        quiet                               = false;
  bool        merge_facets                        = true;
  bool        merge_facets_old_method             = false;
  bool        reject_self_intersections           = false;
  bool        repair_mesh                         = true;
  bool        read_mesh_with_operator             = false;
  bool        debug_input_faces                   = false;
  bool        debug_missing_regions               = false;
  bool        debug_regions                       = false;
  bool        debug_copy_triangulation_into_hole  = false;
  bool        use_new_cavity_algorithm            = true;
  bool        debug_validity                      = false;
  bool        debug_finite_edges_map              = false;
  bool        debug_subconstraints_to_conform     = false;
  bool        debug_verbose_special_cases         = false;
  bool        debug_encroaching_vertices          = false;
  bool        debug_conforming_validation         = false;
  bool        debug_constraint_hierarchy          = false;
  bool        debug_geometric_errors              = false;
  bool        debug_polygon_insertion             = false;
  bool        use_finite_edges_map                = false;
  bool        call_is_valid                       = true;
  double      ratio                               = 0.1;
  double      vertex_vertex_epsilon               = 1e-14;
  double      segment_vertex_epsilon              = 1e-14;
  double      coplanar_polygon_max_angle          = 5.1;
  double      coplanar_polygon_max_distance       = 1e-6;
  std::string failure_assertion_expression        {};
  std::string input_filename                      = CGAL::data_file_path("meshes/mpi.off");
  std::string output_filename                     {"dump.off"};
  std::string dump_patches_after_merge_filename   {};
  std::string dump_surface_mesh_after_merge_filename{};
  std::string dump_patches_borders_prefix         {};
  std::string dump_after_conforming_filename      {};

  CDT_options(int argc, char* argv[]);
};

CDT_options::CDT_options(int argc, char* argv[]) {
  const std::vector<std::string_view> args(argv + 1, argv + argc);
  int positional = 0;

  using std::literals::string_view_literals::operator""sv;
  for (auto it = args.begin(); it != args.end(); ++it) {
    auto get_next_arg_or_error_out = [&it, &args]() -> std::string {
      if(it + 1 == args.end()) {
        error("extra argument required after "sv, *it);
      }
      return std::string(*++it);
    };
    std::string_view arg = *it;
    if(arg == "--merge-facets"sv) {
      merge_facets                        = true;
    } else if(arg == "--no-merge-facets"sv) {
      merge_facets                        = false;
    } else if(arg == "--use-new-cavity-algorithm"sv) {
      use_new_cavity_algorithm            = true;
    } else if(arg == "--use-old-cavity-algorithm"sv) {
      use_new_cavity_algorithm            = false;
    } else if(arg == "--reject-self-intersections"sv) {
      reject_self_intersections           = true;
    } else if(arg == "--no-repair"sv) {
      repair_mesh                         = false;
    } else if(arg == "--read-mesh-with-operator"sv) {
      read_mesh_with_operator             = true;
    } else if(arg == "--merge-facets-old"sv) {
      merge_facets                        = true;
      merge_facets_old_method             = true;
    } else if(arg == "--failure-expression"sv) {
      failure_assertion_expression        = get_next_arg_or_error_out();
    } else if(arg == "--ratio"sv) {
      ratio                               = std::stod(get_next_arg_or_error_out());
    } else if(arg == "--dump-patches-after-merge"sv) {
      dump_patches_after_merge_filename   = get_next_arg_or_error_out();
    } else if(arg == "--dump-patches-borders-prefix"sv) {
      dump_patches_borders_prefix         = get_next_arg_or_error_out();
    } else if(arg == "--dump-surface-mesh-after-merge"sv) {
      dump_surface_mesh_after_merge_filename = get_next_arg_or_error_out();
    } else if(arg == "--dump-after-conforming"sv) {
      dump_after_conforming_filename      = get_next_arg_or_error_out();
    } else if(arg == "--vertex-vertex-epsilon"sv) {
      vertex_vertex_epsilon               = std::stod(get_next_arg_or_error_out());
    } else if(arg == "--segment-vertex-epsilon"sv) {
      segment_vertex_epsilon              = std::stod(get_next_arg_or_error_out());
    } else if(arg == "--coplanar-polygon-max-angle"sv) {
      coplanar_polygon_max_angle          = std::stod(get_next_arg_or_error_out());
    } else if(arg == "--coplanar-polygon-max-distance"sv) {
      coplanar_polygon_max_distance       = std::stod(get_next_arg_or_error_out());
    } else if(arg == "--quiet"sv) {
      quiet                               = true;
    } else if(arg == "--no-is-valid"sv) {
      call_is_valid                       = false;
    } else if(arg == "--debug-input-faces"sv) {
      debug_input_faces                   = true;
    } else if(arg == "--debug-missing-regions"sv) {
      debug_missing_regions               = true;
    } else if(arg == "--debug-regions"sv) {
      debug_regions                       = true;
    } else if(arg == "--debug_copy_triangulation_into_hole"sv) {
      debug_copy_triangulation_into_hole  = true;
    } else if(arg == "--debug-validity"sv) {
      debug_validity                      = true;
    } else if(arg == "--debug-finite-edges-map"sv) {
      debug_finite_edges_map              = true;
    } else if(arg == "--debug-subconstraints-to-conform"sv) {
      debug_subconstraints_to_conform     = true;
    } else if(arg == "--debug-verbose-special-cases"sv) {
      debug_verbose_special_cases         = true;
    } else if(arg == "--debug-encroaching-vertices"sv) {
      debug_encroaching_vertices          = true;
    } else if(arg == "--debug-conforming-validation"sv) {
      debug_conforming_validation         = true;
    } else if(arg == "--debug-constraint-hierarchy"sv) {
      debug_constraint_hierarchy          = true;
    } else if(arg == "--debug-geometric-errors"sv) {
      debug_geometric_errors              = true;
    } else if(arg == "--debug-polygon-insertion"sv) {
      debug_polygon_insertion             = true;
    } else if(arg == "--use-finite-edges-map"sv) {
      use_finite_edges_map                = true;
    } else if(arg == "--verbose"sv || arg == "-V"sv) {
      ++verbose_level;
    } else if(arg == "--help"sv || arg == "-h"sv) {
      need_help                           = true;
    } else if(arg[0] == '-') {
      error("unknown option "sv, arg);
    } else {
      switch(positional) {
        case 0:
          input_filename = arg;
          ++positional;
          break;
        case 1:
          output_filename = arg;
          ++positional;
          break;
        default:
          error("too many arguments"sv);
      }
    }
  }
}

#if NO_TRY_CATCH
# define CDT_3_try      if (true)
# define CDT_3_catch(X) if (false)
# define CDT_3_throw_exception_again
#else
// Else proceed normally.
# define CDT_3_try      try
# define CDT_3_catch(X) catch(X)
# define CDT_3_throw_exception_again throw
#endif

#if CGAL_USE_ITT
#  include <ittnotify.h>
#  define CGAL_CDT_3_TASK_BEGIN(task_handle) \
  std::cerr << "START " << #task_handle << '\n'; \
  __itt_task_begin(cdt_3_domain, __itt_null, __itt_null, task_handle);
#  define CGAL_CDT_3_TASK_END(task_handle) \
  std::cerr << "-STOP " << #task_handle << '\n'; \
  __itt_task_end(cdt_3_domain);

  auto cdt_3_domain = __itt_domain_create("org.cgal.CDT_3");
  auto read_input_task_handle = __itt_string_handle_create("CDT_3: read input file");
  auto merge_facets_task_handle = __itt_string_handle_create("CDT_3: merge facets");
  auto insert_vertices_task_handle = __itt_string_handle_create("CDT_3: insert vertices");
  auto compute_distances_task_handle = __itt_string_handle_create("CDT_3: compute distances");
  auto conforming_task_handle = __itt_string_handle_create("CDT_3: conforming");
  auto cdt_task_handle = __itt_string_handle_create("CDT_3: cdt");
  auto output_task_handle = __itt_string_handle_create("CDT_3: outputs");
  auto validation_task_handle = __itt_string_handle_create("CDT_3: validation");

#else // no ITT
#  define CGAL_CDT_3_TASK_BEGIN(task_handle)
#  define CGAL_CDT_3_TASK_END(task_handle)
#endif // no ITT

void configure_cdt_debug_options(CDT& cdt, const CDT_options& options) {
  cdt.debug_Steiner_points(options.verbose_level > 0);
  cdt.debug_input_faces(options.debug_input_faces);
  cdt.debug_missing_region(options.verbose_level > 1 || options.debug_missing_regions);
  cdt.debug_regions(options.debug_regions);
  cdt.debug_validity(options.debug_validity);
  cdt.debug_finite_edges_map(options.debug_finite_edges_map);
  cdt.debug_subconstraints_to_conform(options.debug_subconstraints_to_conform);
  cdt.debug_verbose_special_cases(options.debug_verbose_special_cases);
  cdt.debug_encroaching_vertices(options.debug_encroaching_vertices);
  cdt.debug_conforming_validation(options.debug_conforming_validation);
  cdt.debug_constraint_hierarchy(options.debug_constraint_hierarchy);
  cdt.debug_geometric_errors(options.debug_geometric_errors);
  cdt.debug_polygon_insertion(options.debug_polygon_insertion);
  cdt.debug_copy_triangulation_into_hole(options.debug_copy_triangulation_into_hole);
  cdt.use_older_cavity_algorithm(!options.use_new_cavity_algorithm);
  cdt.use_finite_edges_map(options.use_finite_edges_map);
  cdt.set_segment_vertex_epsilon(options.segment_vertex_epsilon);
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
    CGAL_CDT_3_TASK_BEGIN(output_task_handle);
    {
      auto dump_tets_to_medit = [](std::string fname,
                              const std::vector<K::Point_3> &points,
                              const std::vector<std::array<std::size_t, 4>> &indexed_tetra,
                              const std::vector<std::size_t> &cell_ids)
      {
        std::ofstream out(fname);
        out.precision(17);
        out << "MeshVersionFormatted 1\nDimension 3\nVertices\n";
        out << points.size() << "\n";
        for (const K::Point_3& p : points)
          out << p << " 0\n";
        out << "Triangles\n0\nTetrahedra\n";
        out << indexed_tetra.size() << "\n";
        for (std::size_t k=0;k<indexed_tetra.size(); ++k)
          out << indexed_tetra[k][0]+1 << " "
                    << indexed_tetra[k][1]+1 << " "
                    << indexed_tetra[k][2]+1 << " "
                    << indexed_tetra[k][3]+1 << " " << cell_ids[k] << "\n";
        out <<"End\n";
      };

      auto& tr = cdt;

      std::unordered_map<CDT::Cell_handle, int /*Subdomain_index*/> cells_map;
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

      std::vector<K::Point_3> points(cdt.number_of_vertices());
      for(auto v: cdt.finite_vertex_handles()) {
        points.at(v->time_stamp() -1) = v->point();
      }
      std::vector<std::array<std::size_t, 4>> indexed_tetra;
      indexed_tetra.reserve(cdt.number_of_cells());
      for(auto ch: cdt.finite_cell_handles()) {
        if(cells_map[ch] > 0) {
          indexed_tetra.push_back({ch->vertex(0)->time_stamp() -1,
                                   ch->vertex(1)->time_stamp() -1,
                                   ch->vertex(2)->time_stamp() -1,
                                   ch->vertex(3)->time_stamp() -1});
        }
      }
      std::vector<std::size_t> cell_idsl(indexed_tetra.size(), 1);
      dump_tets_to_medit(options.output_filename + ".mesh", points, indexed_tetra, cell_idsl);
    }
    {
      std::ofstream dump(options.output_filename);
      dump.precision(17);
#if CGAL_CXX20 && __cpp_lib_concepts >= 201806L && __cpp_lib_ranges >= 201911L
      cdt.write_facets(dump, cdt, std::views::filter(cdt.finite_facets(), [&](auto f) {
          return cdt.is_facet_constrained(f);
      }));
#else
      auto is_facet_constrained = [&](auto f) { return cdt.is_facet_constrained(f); };
      auto it_end = cdt.finite_facets_end();
      cdt.write_facets(dump, cdt,
                        CGAL::make_range(
                          boost::make_filter_iterator(is_facet_constrained,cdt.finite_facets_begin(), it_end),
                          boost::make_filter_iterator(is_facet_constrained,it_end, it_end)));

#endif
    }
    CGAL_CDT_3_TASK_END(output_task_handle);
  };
}

struct Min_distance_result {
  double min_distance;
  std::array<CDT::Vertex_handle, 2> vertices_of_min_edge;
};

Min_distance_result compute_minimum_vertex_distance(const CDT& cdt) {
#if CGAL_CXX20 && __cpp_lib_concepts >= 201806L && __cpp_lib_ranges >= 201911L
  auto [min_sq_distance, min_edge] =
      (std::ranges::min)(cdt.finite_edges() | std::views::transform([&](auto edge) {
                           return std::make_pair(cdt.segment(edge).squared_length(), edge);
                         }));
#else
  auto transform_fct = [&](auto edge) { return std::make_pair(cdt.segment(edge).squared_length(), edge); };
  auto min_p = transform_fct(*cdt.finite_edges_begin());
  for (auto ite=cdt.finite_edges_begin(); ite!=cdt.finite_edges_end(); ++ite)
  {
    auto p = transform_fct(*ite);
    if (p < min_p)
      p = min_p;
  }
  auto [min_sq_distance, min_edge] = min_p;
#endif
  auto min_distance = CGAL::approximate_sqrt(min_sq_distance);
  auto vertices_of_min_edge = cdt.vertices(min_edge);

  return {min_distance, vertices_of_min_edge};
}

void print_minimum_distance_info(const Min_distance_result& min_dist) {
  std::cout << "Min distance between vertices: " << min_dist.min_distance << '\n'
            << "  between vertices:          : "
            << CGAL::IO::oformat(min_dist.vertices_of_min_edge[0], CGAL::With_point_tag{}) << "    "
            << CGAL::IO::oformat(min_dist.vertices_of_min_edge[1], CGAL::With_point_tag{}) << "\n\n";
}

int validate_minimum_vertex_distances(const CDT& cdt, double vertex_vertex_min_distance, const CDT_options& options) {
  auto result = compute_minimum_vertex_distance(cdt);

  if(!options.quiet) {
    print_minimum_distance_info(result);
  }
  if(result.min_distance < vertex_vertex_min_distance) {
    std::cerr << "ERROR: min distance between vertices is too small\n";
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

struct Constraint_distance_result {
  double min_distance;
  CDT::Vertex_handle min_va, min_vb, min_vertex;
};

template<typename BordersOfPatches, typename VertexPointMap>
Constraint_distance_result compute_constraint_vertex_distances_from_patches_borders(
    CDT& cdt,
    const BordersOfPatches& patch_edges,
    const VertexPointMap& tr_vertex_pmap) {

#if CGAL_CXX20 && __cpp_lib_concepts >= 201806L && __cpp_lib_ranges >= 201911L
  auto edge_results = patch_edges
    | std::views::join
    | std::views::transform([&](const auto& edge_pair) {
        auto [vda, vdb] = edge_pair;
        auto va = get(tr_vertex_pmap, vda);
        auto vb = get(tr_vertex_pmap, vdb);
        auto [min_dist, min_v] = cdt.min_distance_and_vertex_between_constraint_and_encroaching_vertex(va, vb);
        return std::make_tuple(CGAL::to_double(min_dist), va, vb, min_v);
      });

  auto min_result = std::ranges::min_element(edge_results, {}, [](const auto& tuple) {
    return std::get<0>(tuple);
  });

  auto [min_distance, min_va, min_vb, min_vertex] = *min_result;
#else
  double min_distance = (std::numeric_limits<double>::max)();
  CDT::Vertex_handle min_va, min_vb, min_vertex;

  for(const auto& edges : patch_edges) {
    for(auto [vda, vdb]: edges) {
      auto va = get(tr_vertex_pmap, vda);
      auto vb = get(tr_vertex_pmap, vdb);
      auto [min_dist, min_v] = cdt.min_distance_and_vertex_between_constraint_and_encroaching_vertex(va, vb);
      if(min_dist < min_distance) {
        min_distance = CGAL::to_double(min_dist);
        min_va = va;
        min_vb = vb;
        min_vertex = min_v;
      }
    }
  }
#endif

  return {min_distance, min_va, min_vb, min_vertex};
}

template<typename VertexPointMap>
Constraint_distance_result compute_constraint_vertex_distances_from_faces(
    CDT& cdt,
    const Mesh& mesh,
    const VertexPointMap& tr_vertex_pmap) {

  double min_distance = (std::numeric_limits<double>::max)();
  CDT::Vertex_handle min_va, min_vb, min_vertex;

  for(auto face_descriptor : faces(mesh)) {
    auto he = halfedge(face_descriptor, mesh);
    const auto end = he;
    do {
      auto va = get(tr_vertex_pmap, source(he, mesh));
      auto vb = get(tr_vertex_pmap, target(he, mesh));
      auto [min_dist, min_v] = cdt.min_distance_and_vertex_between_constraint_and_encroaching_vertex(va, vb);
      if(min_dist < min_distance) {
        min_distance = CGAL::to_double(min_dist);
        min_va = va;
        min_vb = vb;
        min_vertex = min_v;
      }
      he = next(he, mesh);
    } while((he = next(he, mesh)) != end);
  }

  return {min_distance, min_va, min_vb, min_vertex};
}

template<typename BordersOfPatches, typename VertexPointMap>
void validate_constraint_vertex_distances_or_throw(
    CDT& cdt,
    const Mesh& mesh,
    const CDT_options& options,
    const BordersOfPatches& patch_edges,
    const VertexPointMap& tr_vertex_pmap) {

  auto [min_distance, min_va, min_vb, min_vertex] =
      options.merge_facets ? compute_constraint_vertex_distances_from_patches_borders(cdt, patch_edges, tr_vertex_pmap)
                           : compute_constraint_vertex_distances_from_faces(cdt, mesh, tr_vertex_pmap);

  if(!options.quiet) {
    std::cout << "Min distance between constraint segment and vertex: " << min_distance << '\n'
              << "  between segment                                 : "
              << CGAL::IO::oformat(min_va, CDT::Conforming_Dt::with_point) << "    "
              << CGAL::IO::oformat(min_vb, CDT::Conforming_Dt::with_point) << '\n'
              << "  and vertex                                      : "
              << CGAL::IO::oformat(min_vertex, CDT::Conforming_Dt::with_point) << "\n\n";
  }
  cdt.check_segment_vertex_distance_or_throw(min_va, min_vb, min_vertex, min_distance,
                                             CDT::Check_distance::NON_SQUARED_DISTANCE);
}

template<typename PatchIdMap, typename VertexSelectedMap, typename EdgeBorderMap, typename VertexPointMap>
struct Mesh_property_maps {
  PatchIdMap patch_id_map;
  VertexSelectedMap v_selected_map;
  EdgeBorderMap edge_is_border_of_patch_map;
  VertexPointMap mesh_vertex_point_map;
};

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
  int number_of_patches = CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces(
      mesh, pmaps.patch_id_map,
      np::maximum_distance(coplanar_polygon_max_distance)
         .maximum_angle(coplanar_polygon_max_angle));
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
    auto [corner_id_map, corner_id_map_ok] = mesh.add_property_map<vertex_descriptor, std::size_t>("v:corner_id", -1);
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
    double bbox_max_span) {

  int number_of_patches = 0;
  Borders_of_patches patch_edges;

  if(options.merge_facets) {
    CGAL_CDT_3_TASK_BEGIN(merge_facets_task_handle);
    auto start_time = std::chrono::high_resolution_clock::now();

    if(options.merge_facets_old_method) {
      number_of_patches = merge_facets_old_method(mesh, pmaps, number_of_patches);
    } else {
      number_of_patches = merge_facets_region_growing(
          mesh, pmaps, options.coplanar_polygon_max_distance * bbox_max_span,
          options.coplanar_polygon_max_angle, options.dump_surface_mesh_after_merge_filename);
    }
    if (!options.quiet) {
      std::cout << "[timings] detected " << number_of_patches << " patches in " << std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    }
    patch_edges = extract_patch_edges(mesh, pmaps, number_of_patches);
    CGAL_CDT_3_TASK_END(merge_facets_task_handle);
    if(!options.dump_patches_after_merge_filename.empty()) {
      CGAL_CDT_3_TASK_BEGIN(output_task_handle);
      std::ofstream out(options.dump_patches_after_merge_filename);
      CGAL::IO::write_PLY(out, mesh, CGAL::parameters::stream_precision(17));
      CGAL_CDT_3_TASK_END(output_task_handle);
    }
  }

  return patch_edges;
}

template <typename VdToVhPmap, typename BordersOfPatches, typename VertexPointPmap>
void insert_patches_borders_as_constraints(CDT& cdt,
                               BordersOfPatches patch_edges,
                               VdToVhPmap mesh_descriptor_to_vertex_handle_pmap,
                               VertexPointPmap vertex_point_pmap) {

  for(auto& edges : patch_edges) {
    if(edges.empty())
      continue;
    auto polylines = CGAL::segment_soup_to_polylines(edges);
    while(true) {
      const auto non_closed_polylines_begin =
          std::partition(polylines.begin(), polylines.end(),
                         [](const auto& polyline) { return polyline.front() == polyline.back(); });
      if(non_closed_polylines_begin == polylines.end())
        break;
      edges.clear();
      for(auto it = non_closed_polylines_begin; it != polylines.end(); ++it) {
        auto& polyline = *it;
        for(auto it = polyline.begin(), end = polyline.end() - 1; it != end; ++it) {
          edges.emplace_back(*it, *(it + 1));
        }
      }
      polylines.erase(non_closed_polylines_begin, polylines.end());
      auto other_polylines = CGAL::segment_soup_to_polylines(edges);
      polylines.insert(polylines.end(),
                       std::make_move_iterator(other_polylines.begin()),
                       std::make_move_iterator(other_polylines.end()));
    }

    if(polylines.size() > 1) {
      double max_sq_length = 0;
      auto longest_it = polylines.begin();
      for(auto it = polylines.begin(); it != polylines.end(); ++it) {
        auto& polyline = *it;
        using CGAL::Bbox_3;
        Bbox_3 bb;
        for(auto v : polyline) {
          bb = bb + Bbox_3(get(vertex_point_pmap, v).bbox());
        }
        double sq_diagonal_length = CGAL::square(bb.xmax() - bb.xmin()) +
                                    CGAL::square(bb.ymax() - bb.ymin()) +
                                    CGAL::square(bb.zmax() - bb.zmin());
        if(sq_diagonal_length > max_sq_length) {
          max_sq_length = sq_diagonal_length;
          longest_it = it;
        }
      }
      if(longest_it != polylines.begin()) {
        std::iter_swap(longest_it, polylines.begin());
      }
    }

    std::optional<int> face_index;
    for(auto& polyline : polylines) {
      CGAL_assertion(!polyline.empty() && polyline.front() == polyline.back());
      polyline.pop_back();
      auto range_of_vertices = CGAL::make_transform_range_from_property_map(polyline, mesh_descriptor_to_vertex_handle_pmap);
      face_index = cdt.insert_constrained_face(range_of_vertices, false,
                                               face_index ? *face_index : -1);
    }
  }
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
  CDT cdt;
  configure_cdt_debug_options(cdt, options);

  auto bbox = compute_bounding_box(mesh, options);
  auto bbox_max_span = (std::max)(bbox.x_span(), (std::max)(bbox.y_span(), bbox.z_span()));

  auto pmaps = setup_mesh_property_maps(mesh);
  auto patch_edges = maybe_merge_facets(mesh, options, pmaps, bbox_max_span);
  if(!options.dump_patches_borders_prefix.empty()) {
    CGAL_CDT_3_TASK_BEGIN(output_task_handle);
    dump_patches_borders(patch_edges, pmaps, options.dump_patches_borders_prefix);
    CGAL_CDT_3_TASK_END(output_task_handle);
  }

  int exit_code = EXIT_SUCCESS;

  auto [mesh_descriptor_to_vertex_handle_pmap, tr_vertex_pmap_ok] = mesh.add_property_map<vertex_descriptor, CDT::Vertex_handle>("tr_vertex");
  assert(tr_vertex_pmap_ok); CGAL_USE(tr_vertex_pmap_ok);

  CGAL_CDT_3_TASK_BEGIN(insert_vertices_task_handle);
  auto start_time = std::chrono::high_resolution_clock::now();
  CDT::Cell_handle hint{};
  for(auto v: vertices(mesh)) {
    if(options.merge_facets && false == get(pmaps.v_selected_map, v)) continue;
    auto vh = cdt.insert(get(pmaps.mesh_vertex_point_map, v), hint, false);
    hint = vh->cell();
    put(mesh_descriptor_to_vertex_handle_pmap, v, vh);
  }
  if(!options.quiet) {
    std::cout << "[timings] inserted vertices in " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    std::cout << "Number of vertices: " << cdt.number_of_vertices() << "\n\n";
  }
  if(cdt.dimension() < 3) {
    if(!options.quiet) {
      std::cout << "current is 2D... inserting the 8 vertices of an extended bounding box\n";
    }

    auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
    auto dx = bbox.x_span();
    auto dy = bbox.y_span();
    auto dz = bbox.z_span();
    auto max_span = (std::max)(dx, (std::max)(dy, dz));
    if(dx == 0) dx = max_span;
    if(dy == 0) dy = max_span;
    if(dz == 0) dz = max_span;

    cdt.insert(Point(bbox.xmin() - dx, bbox.ymin() - dy, bbox.zmin() - dz));
    cdt.insert(Point(bbox.xmin() - dx, bbox.ymax() + dy, bbox.zmin() - dz));
    cdt.insert(Point(bbox.xmin() - dx, bbox.ymin() - dy, bbox.zmax() + dz));
    cdt.insert(Point(bbox.xmin() - dx, bbox.ymax() + dy, bbox.zmax() + dz));
    cdt.insert(Point(bbox.xmax() + dx, bbox.ymin() - dy, bbox.zmin() - dz));
    cdt.insert(Point(bbox.xmax() + dx, bbox.ymax() + dy, bbox.zmin() - dz));
    cdt.insert(Point(bbox.xmax() + dx, bbox.ymin() - dy, bbox.zmax() + dz));
    cdt.insert(Point(bbox.xmax() + dx, bbox.ymax() + dy, bbox.zmax() + dz));
  }
  CGAL_CDT_3_TASK_END(insert_vertices_task_handle);

  start_time = std::chrono::high_resolution_clock::now();
  CGAL_CDT_3_TASK_BEGIN(compute_distances_task_handle);

  exit_code = validate_minimum_vertex_distances(cdt, options.vertex_vertex_epsilon * bbox_max_span, options);
  if(exit_code != EXIT_SUCCESS) {
    return exit_code;
  }

  validate_constraint_vertex_distances_or_throw(cdt, mesh, options, patch_edges, mesh_descriptor_to_vertex_handle_pmap);
  CGAL_CDT_3_TASK_END(compute_distances_task_handle);

  if(!options.quiet) {
    std::cout << "[timings] compute distances on " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
  }

  CGAL_CDT_3_TASK_BEGIN(conforming_task_handle);

  int poly_id = 0;
  auto output_on_exit_scope_guard = CGAL::make_scope_exit(create_output_finalizer(cdt, options));
  start_time = std::chrono::high_resolution_clock::now();
  if(options.merge_facets) {
    insert_patches_borders_as_constraints(cdt, std::move(patch_edges), mesh_descriptor_to_vertex_handle_pmap,
                                          pmaps.mesh_vertex_point_map);
  } else {
    for(auto face_descriptor : faces(mesh)) {
      std::vector<Point_3> polygon;
      const auto he = halfedge(face_descriptor, mesh);
      for(auto vertex_it : CGAL::vertices_around_face(he, mesh)) {
        polygon.push_back(get(pmaps.mesh_vertex_point_map, vertex_it));
      }
      if(cdt.debug_polygon_insertion()) {
        std::cerr << "NEW POLYGON #" << poly_id << '\n';
      }
      try {
        [[maybe_unused]] auto id = cdt.insert_constrained_polygon(polygon, false);
        assert(id == poly_id);
        ++poly_id;
      } catch(int error) {
        exit_code = error;
      }
      // std::ofstream dump("dump.binary.cgal");
      // CGAL::Mesh_3::save_binary_file(dump, cdt);
    }
  } // not merge_facets
  if(!options.quiet) {
    std::cout << "[timings] registered facets in " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
  }
  start_time = std::chrono::high_resolution_clock::now();
  cdt.restore_Delaunay();
  if(!options.quiet) {
    std::cout << "[timings] restored Delaunay (conforming of facets borders) in " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
  }

  CGAL_CDT_3_TASK_END(conforming_task_handle);

  if(!options.dump_after_conforming_filename.empty()) {
    CGAL_CDT_3_TASK_BEGIN(output_task_handle);
    using Vertex_index = Mesh::Vertex_index;
    [[maybe_unused]] std::size_t time_stamp_counter = 0u;
    for(auto v: cdt.finite_vertex_handles()) {
      [[maybe_unused]] const auto time_stamp = v->time_stamp();
      assert(++time_stamp_counter == time_stamp);
      if(!v->ccdt_3_data().is_Steiner_vertex_on_edge()) continue;
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
    CGAL_CDT_3_TASK_END(output_task_handle);
  }

  if(!options.quiet) {
    std::cout << "Number of vertices after conforming: " << cdt.number_of_vertices() << "\n\n";
  }
  CGAL_CDT_3_TASK_BEGIN(validation_task_handle);
  CGAL_assertion(!options.call_is_valid || cdt.Base_triantulation::is_valid(true));
  CGAL_assertion(!options.call_is_valid || cdt.is_valid(true));
  CGAL_assertion(!options.call_is_valid || cdt.is_conforming());
  CGAL_CDT_3_TASK_END(validation_task_handle);
  if(exit_code == EXIT_SUCCESS) {
    try {
      CGAL_CDT_3_TASK_BEGIN(cdt_task_handle);
      start_time = std::chrono::high_resolution_clock::now();
      cdt.restore_constrained_Delaunay();
      if(!options.quiet) {
        std::cout << "[timings] restored constrained Delaunay in " << std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
        std::cout << "Number of vertices after CDT: " << cdt.number_of_vertices() << "\n\n";
      }
      CGAL_CDT_3_TASK_END(cdt_task_handle);
    } catch(int error) {
      exit_code = error;
    }
  }

  CGAL_CDT_3_TASK_BEGIN(validation_task_handle);
  CGAL_assertion(!options.call_is_valid || cdt.is_conforming());
  CGAL_assertion(!options.call_is_valid || cdt.is_valid(true));
  CGAL_CDT_3_TASK_END(validation_task_handle);

  return exit_code;
}

int bisect_errors(Mesh mesh, CDT_options options) {
  auto nb_buckets = static_cast<int>(std::floor(1 / options.ratio)) + 1;
  std::cerr << "RATIO: " << options.ratio << '\n';

  const Mesh orig_mesh{mesh};
  Mesh bad_mesh{mesh};

  int exit_code = EXIT_SUCCESS;

  for(int bucket = 0; bucket < nb_buckets;) {
    const auto nb_faces = mesh.number_of_faces();
    auto nb_to_skip = static_cast<int>(std::round(nb_faces * options.ratio));
    if(nb_to_skip < 1) {
      nb_to_skip = 1;
      nb_buckets = nb_faces;
    }
    if(bucket == 0) {
      std::cerr << "NB BUCKETS: " << nb_buckets << '\n';
    }

    auto simplify = [&](Mesh& m) {
      std::cerr << "nb_to_skip: " << nb_to_skip << '\n';
      std::cerr << "bucket: " << bucket << '\n';
      const auto start = (std::min)(bucket * nb_to_skip, static_cast<int>(m.number_of_faces()));
      const auto end = (std::min)(start + nb_to_skip, static_cast<int>(m.number_of_faces()));
      std::cerr << "SKIP from " << start << " to " << end << '\n';
      for(auto i = end - 1; i >= start; --i) {
        const auto f = m.faces().begin() + i;
        CGAL::Euler::remove_face(halfedge(*f, m), m);
      }
      assert(m.is_valid(true));
      std::cerr << "number of faces: " << m.number_of_faces() << '\n';
      if(m.number_of_faces() >= nb_faces) {
        std::cerr << "ERROR: could not simplify mesh\n";
        return false;
      }
      return true;
    };

    if(simplify(mesh)) {
      std::ofstream current("current_mesh.off");
      current.precision(17);
      current << mesh;
      current.close();

      try {
        auto code = go(mesh, options);
        if(code != EXIT_SUCCESS) {
          exit_code = code;
        }
      } catch(std::exception& e) {
        std::cerr << "CAUGHT EXCEPTION: " << e.what() << '\n';
        if(std::string(e.what()).find(options.failure_assertion_expression) != std::string::npos) {
          exit_code = EXIT_FAILURE;
          std::cerr << "BAD MESH! " << mesh.number_of_faces() << " faces\n";
          std::ofstream bad("bad_mesh.off");
          bad.precision(17);
          bad << mesh;
          bad_mesh = mesh;
          bucket = 0;
          continue;
        } else {
          exit_code = EXIT_FAILURE;
          std::cerr << "ERROR MESH: " << e.what() << '\n';
          std::ofstream error("error_mesh.off");
          error.precision(17);
          error << mesh;
          std::cerr << "go on...\n";
        }
      }
      std::cerr << "GOOD MESH :-( " << mesh.number_of_faces() << " faces\n";
    }
    mesh = bad_mesh;
    ++bucket;
  }
  if(bad_mesh.number_of_faces() < orig_mesh.number_of_faces()) {
    std::cerr << "FINAL BAD MESH: " << bad_mesh.number_of_faces() << " faces\n";
    std::ofstream final_bad("final_bad_mesh.off");
    final_bad.precision(17);
    final_bad << bad_mesh;
    return go(bad_mesh, options);
  }
  return exit_code;
}

int main(int argc, char* argv[]) {
  CDT::Conforming_Dt::with_offset.offset = -1;
  CDT::Conforming_Dt::with_point.offset = -1;
  CDT::Conforming_Dt::with_point_and_info.offset = -1;
  std::cerr.precision(17);
  std::cout.precision(17);

  CDT_options options(argc, argv);
  if(options.need_help) {
    help(std::cout);
    return 0;
  }

  CGAL_CDT_3_TASK_BEGIN(read_input_task_handle);
  auto start_time = std::chrono::high_resolution_clock::now();

  CGAL::CDT_3_read_polygon_mesh_output<Mesh> read_mesh_result;
  if(options.read_mesh_with_operator) {
    std::ifstream in(options.input_filename);
    if(!in) {
      std::cerr << "Cannot open file: " << options.input_filename << std::endl;
      return EXIT_FAILURE;
    }
    Mesh mesh;
    in >> mesh;
    if(!in) {
      std::cerr << "Error reading mesh with operator>>" << std::endl;
      return EXIT_FAILURE;
    }
    read_mesh_result.polygon_mesh = std::move(mesh);
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
  Mesh mesh = std::move(*read_mesh_result.polygon_mesh);
  if(!options.quiet) {
    std::cout << "[timings] read mesh in " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    std::cout << "Number of vertices: " << mesh.number_of_vertices() << '\n';
    std::cout << "Number of edges: " << mesh.number_of_edges() << '\n';
    std::cout << "Number of faces: " << mesh.number_of_faces() << "\n\n";

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
  CGAL_CDT_3_TASK_END(read_input_task_handle);

  if(options.reject_self_intersections && read_mesh_result.polygon_soup_self_intersects) {
    std::cerr << "ERROR: input mesh self-intersects\n";
    return EXIT_FAILURE;
  }

  if(!options.failure_assertion_expression.empty()) {
    return bisect_errors(std::move(mesh), options);
  }

  auto exit_code = go(std::move(mesh), options);
  if(!options.quiet) {
    std::cout << "[timings] total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    if(exit_code != 0) std::cout << "ERROR with exit code " << exit_code << '\n';
  }
  return exit_code;
}
