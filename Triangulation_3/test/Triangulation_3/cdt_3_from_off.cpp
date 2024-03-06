//#define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#define NO_TRY_CATCH 1
// #define CGAL_DEBUG_CDT_3 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/region_growing.h>

#include <boost/graph/graph_traits.hpp>

#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <ranges>
#include <optional>
#include <chrono>

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

#if CGAL_CDT_3_USE_EPECK

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel;

#else // use Epick

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

#endif // use Epick

using Vb = CGAL::Constrained_Delaunay_triangulation_vertex_base_3<K>;
using Cb = CGAL::Constrained_Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
using CDT = CGAL::Constrained_Delaunay_triangulation_3<Delaunay>;
using Point = Delaunay::Point;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using edge_descriptor   = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor   = boost::graph_traits<Mesh>::face_descriptor;

struct CDT_options
{
  int verbose = 0;
  bool quiet = false;
  bool merge_facets = true;
  bool merge_facets_old_method = false;
  bool repair_mesh = true;
  bool debug_missing_regions = false;
  bool debug_regions = false;
  bool debug_copy_triangulation_into_hole = false;
  bool use_new_cavity_algorithm = true;
  bool debug_validity = false;
  double ratio = 0.1;
  double vertex_vertex_epsilon = 1e-14;
  double segment_vertex_epsilon = 1e-14;
  std::string failure_assertion_expression{};
  std::string input_filename = CGAL::data_file_path("meshes/mpi.off");
  std::string output_filename{"dump.off"};
  std::string dump_patches_after_merge_filename{};
  std::string dump_patches_borders_prefix{};
  std::string dump_after_conforming_filename{};
};

int go(Mesh, CDT_options);

void help(std::ostream& out) {
  out << R"(
Usage: cdt_3_from_off [options] input.off output.off

  input.off: input mesh
  output.off: output mesh

  --merge-facets/--no-merge-facets: merge facets into patches (set by default)
  --use-new-cavity-algorithm/--use-old-cavity-algorithm: use new or old cavity algorithm (default: new)
  --failure-expression <expression>: expression to detect bad meshratio)
  --ratio <double>: ratio of faces to remove (default: 0.1), if --failure-expression is used
  --vertex-vertex-epsilon <double>: epsilon for vertex-vertex min distance (default: 1e-6)
  --segment-vertex-epsilon <double>: epsilon for segment-vertex min distance (default: 0)

  --dump-patches-after-merge <filename.ply>: dump patches after merging facets in PLY
  --dump-patches-borders-prefix <filenames_prefix>: dump patches borders
  --dump-after-conforming <filename.off>: dump mesh after conforming in OFF

  -V/--verbose: verbose (can be used several times)
  --debug-missing-regions: debug missing regions
  --debug-regions: debug regions
  --debug_copy_triangulation_into_hole: debug copy_triangulation_into_hole
  --debug-validity: add is_valid checks after modifications to the TDS

  --quiet: do not print anything
  --help: print this help
)";
}

int main(int argc, char* argv[])
{
  CDT::Conforming_Dt::with_offset.offset = -1;
  CDT::Conforming_Dt::with_point.offset = -1;
  CDT::Conforming_Dt::with_point_and_info.offset = -1;
  std::cerr.precision(17);
  std::cout.precision(17);

  CDT_options options;
  int positional = 0;

  for(int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if(arg == "--merge-facets") {
      options.merge_facets = true;
    } else if(arg == "--no-merge-facets") {
      options.merge_facets = false;
    } else if(arg == "--use-new-cavity-algorithm") {
      options.use_new_cavity_algorithm = true;
    } else if(arg == "--use-old-cavity-algorithm") {
      options.use_new_cavity_algorithm = false;
    } else if(arg == "--no-repair") {
      options.repair_mesh = false;
    } else if(arg == "--merge-facets-old") {
      options.merge_facets = true;
      options.merge_facets_old_method = true;
    } else if(arg == "--failure-expression") {
      assert(i + 1 < argc);
      options.failure_assertion_expression = argv[++i];
    } else if(arg == "--ratio") {
      assert(i + 1 < argc);
      options.ratio = std::stod(argv[++i]);
    } else if(arg == "--dump-patches-after-merge") {
      assert(i + 1 < argc);
      options.dump_patches_after_merge_filename = argv[++i];
    } else if(arg == "--dump-patches-borders-prefix") {
      assert(i + 1 < argc);
      options.dump_patches_borders_prefix = argv[++i];
    } else if(arg == "--dump-after-conforming") {
      assert(i + 1 < argc);
      options.dump_after_conforming_filename = argv[++i];
    } else if(arg == "--vertex-vertex-epsilon") {
      assert(i + 1 < argc);
      options.vertex_vertex_epsilon = std::stod(argv[++i]);
    } else if(arg == "--segment-vertex-epsilon") {
      assert(i + 1 < argc);
      options.segment_vertex_epsilon = std::stod(argv[++i]);
    } else if(arg == "--quiet") {
      options.quiet = true;
    } else if(arg == "--debug-missing-regions") {
      options.debug_missing_regions = true;
    } else if(arg == "--debug-regions") {
      options.debug_regions = true;
    } else if(arg == "--debug_copy_triangulation_into_hole") {
      options.debug_copy_triangulation_into_hole = true;
    } else if(arg == "--debug-validity") {
      options.debug_validity = true;
    } else if(arg == "-V") {
      ++options.verbose;
    } else if(arg == "--help") {
      help(std::cout);
      return 0;
    } else if(arg[0] == '-') {
      std::cerr << "Unknown option: " << arg << '\n';
      help(std::cerr);
      return 1;
    } else {
      switch(positional) {
        case 0:
          options.input_filename = arg;
          ++positional;
          break;
        case 1:
          options.output_filename = arg;
          ++positional;
          break;
        case 2:
          options.ratio = std::stod(arg);
          ++positional;
          break;
        default:
          std::cerr << "Too many arguments\n";
          return 1;
      }
    }
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  Mesh mesh;
  const bool ok = options.repair_mesh
                      ? CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(options.input_filename, mesh)
                      : CGAL::IO::read_polygon_mesh(options.input_filename, mesh);
  if(!ok) {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }
  if(!options.quiet) {
    std::cout << "[timings] read mesh in " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    std::cout << "Number of vertices: " << mesh.number_of_vertices() << '\n';
    std::cout << "Number of edges: " << mesh.number_of_edges() << '\n';
    std::cout << "Number of faces: " << mesh.number_of_faces() << "\n\n";
  }

  if(options.failure_assertion_expression.empty()) {
    auto exit_code = go(std::move(mesh), std::move(options));
    if(!options.quiet) {
      std::cout << "[timings] total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    }
    return exit_code;
  }
  auto nb_buckets = static_cast<int>(std::floor(1 / options.ratio)) + 1;
  std::cerr << "RATIO: " << options.ratio << '\n';

  const Mesh orig_mesh{mesh};
  Mesh bad_mesh{mesh};
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
        std::exit(EXIT_FAILURE);
      }
    };

    simplify(mesh);
    std::ofstream current("current_mesh.off");
    current.precision(17);
    current << mesh;
    current.close();

    try {
      go(mesh, options);
    } catch(CGAL::Failure_exception& e) {
      std::cerr << "CAUGHT EXCEPTION: " << e.what() << '\n';
      if(std::string(e.what()).find(options.failure_assertion_expression) != std::string::npos)
      {
        std::cerr << "BAD MESH! " << mesh.number_of_faces() << " faces\n";
        std::ofstream bad("bad_mesh.off");
        bad.precision(17);
        bad << mesh;
        bad_mesh = mesh;
        bucket = 0;
        continue;
      } else {
        std::cerr << "ERROR MESH: " << e.what() << '\n';
        std::ofstream error("error_mesh.off");
        error.precision(17);
        error << mesh;
        std::cerr << "go on...\n";
      }
    }
    std::cerr << "GOOD MESH :-( " << mesh.number_of_faces() << " faces\n";
    mesh = bad_mesh;
    ++bucket;
  }
}

template <typename Range_of_segments>
auto segment_soup_to_polylines(Range_of_segments&& segment_soup) {
  using Point = decltype([&](){ using std::begin; auto [a, b] = *begin(segment_soup); return a; } ());

  std::vector<std::vector<Point>> polylines;

  using Graph = boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Point>;
  using Map2v = std::map<Point, typename Graph::vertex_descriptor>;
  Graph graph;
  Map2v map2v;
  auto get_v = [&](const Point& p) {
    auto it = map2v.find(p);
    if(it != map2v.end()) return it->second;
    auto v = boost::add_vertex(p, graph);
    map2v.emplace(p, v);
    return v;
  };
  for(auto [a, b]: segment_soup) {
    auto va = get_v(a);
    auto vb = get_v(b);
    boost::add_edge(va, vb, graph);
  }

  struct Polylines_visitor
  {
    Graph& graph;
    std::vector<std::vector<Point>>& polylines;

    void start_new_polyline() { polylines.emplace_back(); }
    void add_node(typename Graph::vertex_descriptor vd) { polylines.back().push_back(graph[vd]); }
    void end_polyline() {}
  };
  Polylines_visitor visitor{graph, polylines};
  CGAL::split_graph_into_polylines(graph, visitor);

  return polylines;
}

int go(Mesh mesh, CDT_options options) {
  CDT cdt;
  cdt.debug_Steiner_points(options.verbose > 0);
  cdt.debug_missing_region(options.verbose > 1 || options.debug_missing_regions);
  cdt.debug_regions(options.debug_regions);
  cdt.debug_validity(options.debug_validity);
  cdt.debug_copy_triangulation_into_hole(options.debug_copy_triangulation_into_hole);
  cdt.use_older_cavity_algorithm(!options.use_new_cavity_algorithm);
  cdt.set_segment_vertex_epsilon(options.segment_vertex_epsilon);

  const auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  double d_x = bbox.xmax() - bbox.xmin();
  double d_y = bbox.ymax() - bbox.ymin();
  double d_z = bbox.zmax() - bbox.zmin();

  const double bbox_max_width = (std::max)(d_x, (std::max)(d_y, d_z));

  double epsilon = options.vertex_vertex_epsilon;

  if(!options.quiet) {
    std::cout << "Bbox width                   : " << bbox_max_width << '\n'
              << "Epsilon                      : " << epsilon << '\n'
              << "Epsilon * Bbox width         : " << epsilon * bbox_max_width << "\n\n";
  }

  auto pmap = get(CGAL::vertex_point, mesh);

  auto [patch_id_map, patch_id_map_ok] = mesh.add_property_map<face_descriptor, int>("f:patch_id", -2);
  assert(patch_id_map_ok); CGAL_USE(patch_id_map_ok);
  auto [v_selected_map, v_selected_map_ok] = mesh.add_property_map<vertex_descriptor, bool>("v:selected", false);
  assert(v_selected_map_ok); CGAL_USE(v_selected_map_ok);
  int nb_patches = 0;
  std::vector<std::vector<std::pair<vertex_descriptor, vertex_descriptor>>> patch_edges;
  if(options.merge_facets) {
    auto start_time = std::chrono::high_resolution_clock::now();

    if(options.merge_facets_old_method) {
      for(auto f: faces(mesh))
      {
        if(get(patch_id_map, f) >= 0) continue;
        std::stack<face_descriptor> f_stack;
        f_stack.push(f);
        while(!f_stack.empty()) {
          auto f = f_stack.top();
          f_stack.pop();
          if(get(patch_id_map, f) >= 0) continue;
          put(patch_id_map, f, nb_patches);
          for(auto h: CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
            auto opp = opposite(h, mesh);
            if(is_border_edge(opp, mesh)) {
              continue;
            }
            auto n = face(opp, mesh);
            auto a = get(pmap, source(h, mesh));
            auto b = get(pmap, target(h, mesh));
            auto c = get(pmap, target(next(h, mesh), mesh));
            auto d = get(pmap, target(next(opp, mesh), mesh));
            if(CGAL::orientation(a, b, c, d) != CGAL::COPLANAR) {
              continue;
            }
            if(get(patch_id_map, n) >= 0) continue;
            f_stack.push(n);
          }
        }
        ++nb_patches;
      }
    } else {
      namespace np = CGAL::parameters;
      nb_patches = CGAL::Polygon_mesh_processing::region_growing_of_planes_on_faces(
          mesh, patch_id_map, np::maximum_distance(epsilon * bbox_max_width).maximum_angle(1));
    }
    for(auto f: faces(mesh)) {
      if(get(patch_id_map, f) < 0) {
        std::cerr << "warning: face " << f << " has no patch id! Reassign it to " << nb_patches << '\n';
        for(auto h: CGAL::halfedges_around_face(halfedge(f, mesh), mesh)) {
          std::cerr << "  " << target(h, mesh) << ", point " << mesh.point(target(h, mesh)) << '\n';
        }
        put(patch_id_map, f, nb_patches++);
      }
    }
    if (!options.quiet) {
      std::cout << "[timings] detected " << nb_patches << " patches in " << std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
    }
    patch_edges.resize(nb_patches);
    for(auto h: halfedges(mesh))
    {
      if(is_border(h, mesh)) continue;
      auto f = face(h, mesh);
      auto patch_id = get(patch_id_map, f);
      auto opp = opposite(h, mesh);
      if(is_border(opp, mesh) || patch_id != get(patch_id_map, face(opp, mesh))) {
        auto va = source(h, mesh);
        auto vb = target(h, mesh);
        patch_edges[patch_id].emplace_back(va, vb);
        put(v_selected_map, va, true);
        put(v_selected_map, vb, true);
      }
    }
    if(!options.dump_patches_after_merge_filename.empty()) {
      std::ofstream out(options.dump_patches_after_merge_filename);
      CGAL::IO::write_PLY(out, mesh);
    }
  }
  if(!options.dump_patches_borders_prefix.empty()) {
    for(int i = 0; i < nb_patches; ++i) {
      std::stringstream ss;
      ss << options.dump_patches_borders_prefix << i << ".polylines.txt";
      std::ofstream out(ss.str());
      out.precision(17);
      const auto& edges = patch_edges[i];
      std::cerr << "Patch p#" << i << " has " << edges.size() << " edges\n";
      const auto polylines = segment_soup_to_polylines(edges);
      for(const auto& polyline: polylines) {
        out << polyline.size() << "    ";
        for(auto v: polyline) {
          out << get(pmap, v) << "  ";
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
  }

  int exit_code = EXIT_SUCCESS;

  auto finally = [&cdt, &options]() {
    {
      std::ofstream dump("dump.binary.cgal");
      CGAL::IO::save_binary_file(dump, cdt);
    }
    {
      std::ofstream dump(options.output_filename);
      dump.precision(17);
      cdt.write_facets(dump, cdt, std::views::filter(cdt.finite_facets(), [&](auto f) {
          return cdt.is_constrained(f);
      }));
    }
    {
      std::ofstream missing_faces("dump_missing_faces.polylines.txt");
      missing_faces.precision(17);
      cdt.recheck_constrained_Delaunay();
      if(cdt.write_missing_subfaces_file(missing_faces)) {
        std::cerr << "ERROR: Missing subfaces!\n";
      }
    }
    {
      std::ofstream missing_edges("dump_missing_segments.polylines.txt");
      missing_edges.precision(17);
      if(cdt.write_missing_segments_file(missing_edges)) {
        std::cerr << "ERROR: Missing segments!\n";
      }
    }
  };

  auto [tr_vertex_pmap, tr_vertex_pmap_ok] = mesh.add_property_map<vertex_descriptor, CDT::Vertex_handle>("tr_vertex");
  assert(tr_vertex_pmap_ok); CGAL_USE(tr_vertex_pmap_ok);

  auto start_time = std::chrono::high_resolution_clock::now();
  for(auto v: vertices(mesh)) {
    if(options.merge_facets && false == get(v_selected_map, v)) continue;
    auto vh = cdt.insert(get(pmap, v));
    put(tr_vertex_pmap, v, vh);
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
    if(d_x == 0) d_x = bbox_max_width;
    if(d_y == 0) d_y = bbox_max_width;
    if(d_z == 0) d_z = bbox_max_width;

    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
    cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));
  }
  start_time = std::chrono::high_resolution_clock::now();
  {
    auto [min_sq_distance, min_edge] = std::ranges::min(
        cdt.finite_edges() | std::views::transform([&](auto edge) { return std::make_pair(cdt.segment(edge).squared_length(), edge); }));
    auto min_distance = CGAL::approximate_sqrt(min_sq_distance);
    auto vertices_of_min_edge = cdt.vertices(min_edge);
    if(!options.quiet) {
      std::cout << "Min distance between vertices: " << min_distance << '\n'
                << "  between vertices:          : " << CGAL::IO::oformat(vertices_of_min_edge[0], CGAL::With_point_tag{})
                << "    " << CGAL::IO::oformat(vertices_of_min_edge[1], CGAL::With_point_tag{}) << "\n\n";
    }
    if(min_distance < epsilon * bbox_max_width) {
      std::cerr << "ERROR: min distance between vertices is too small\n";
      exit_code = EXIT_FAILURE;
      return exit_code;
    }
  }
  {
    double min_distance = std::numeric_limits<double>::max();
    CDT::Vertex_handle min_va, min_vb, min_vertex;
    if(options.merge_facets) {
      for(int i = 0; i < nb_patches; ++i) {
        const auto& edges = patch_edges[i];
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
    } else {
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
    }
    if(!options.quiet) {
      std::cout << "Min distance between constraint segment and vertex: " << min_distance << '\n'
                << "  between segment                                 : "
                << CGAL::IO::oformat(min_va, CDT::Conforming_Dt::with_point) << "    "
                << CGAL::IO::oformat(min_vb, CDT::Conforming_Dt::with_point) << '\n'
                << "  and vertex:                                     : "
                << CGAL::IO::oformat(min_vertex, CDT::Conforming_Dt::with_point) << "\n\n";
    }
    cdt.check_segment_vertex_distance_or_throw(min_va, min_vb, min_vertex, min_distance,
                                               CDT::Check_distance::NON_SQUARED_DISTANCE);
  }
  if(!options.quiet) {
    std::cout << "[timings] compute distances on " << std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
  }
  int poly_id = 0;
  CDT_3_try {
    start_time = std::chrono::high_resolution_clock::now();
    if(options.merge_facets) {
      for(int i = 0; i < nb_patches; ++i) {
        auto& edges = patch_edges[i];
        auto polylines = segment_soup_to_polylines(edges);
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
          auto other_polylines = segment_soup_to_polylines(edges);
          polylines.insert(polylines.end(),
                           std::make_move_iterator(other_polylines.begin()),
                           std::make_move_iterator(other_polylines.end()));
        }

        std::optional<int> face_index;
        for(auto& polyline: polylines) {
          assert(polyline.front() == polyline.back());
          polyline.pop_back();
          if(face_index) {
            cdt.insert_constrained_polygon(
              polyline | std::views::transform([&](vertex_descriptor v) { return get(pmap, v); }),
              false,
              *face_index);
          } else {
            face_index = cdt.insert_constrained_polygon(
              polyline | std::views::transform([&](vertex_descriptor v) { return get(pmap, v); }),
              false);
          }
        }
      }
    } else {
      for(auto face_descriptor : faces(mesh)) {
        std::vector<Point_3> polygon;
        const auto he = halfedge(face_descriptor, mesh);
        for(auto vertex_it : CGAL::vertices_around_face(he, mesh)) {
          polygon.push_back(get(pmap, vertex_it));
        }
  #if CGAL_DEBUG_CDT_3
        std::cerr << "NEW POLYGON #" << poly_id << '\n';
  #endif // CGAL_DEBUG_CDT_3
        const auto coplanar = polygon.size() < 3 ||
            std::all_of(polygon.begin(), polygon.end(),
                        [p1 = polygon[0], p2 = polygon[1], p3 = polygon[2]](auto p) {
                          const auto coplanar =
                              CGAL::orientation(p1, p2, p3, p) == CGAL::COPLANAR;
                          if(!coplanar) {
                            std::cerr << "Non coplanar points: " << p1 << ", " << p2
                                      << ", " << p3 << ", " << p << '\n'
                                      << "  volume: " << volume(p1, p2, p3, p) << '\n';

                          }
                          return coplanar;
                        });
        if(!coplanar) {
          std::ofstream out(std::string("dump_noncoplanar_polygon_") + std::to_string(poly_id) + ".off");
          out.precision(17);
          out << "OFF\n" << polygon.size() << " 1 0\n";
          for(auto p : polygon) {
            out << p << '\n';
          }
          out << polygon.size() << ' ';
          for(std::size_t i = 0u, end = polygon.size(); i < end; ++i) {
            out << ' ' << i;
          }
          out << '\n';
          std::cerr << "Polygon is not coplanar\n";
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

    if(!options.dump_after_conforming_filename.empty()) {
      using Vertex_index = Mesh::Vertex_index;
      [[maybe_unused]] std::size_t time_stamp_counter = 0u;
      for(auto v: cdt.finite_vertex_handles()) {
        const auto time_stamp = v->time_stamp();
        assert(++time_stamp_counter == time_stamp);
        if(!v->is_Steiner_vertex_on_edge()) continue;
        const auto [va, vb] = cdt.ancestors_of_Steiner_vertex_on_edge(v);
        const auto index_va = Vertex_index{static_cast<unsigned>(va->time_stamp() - 1)};
        const auto index_vb = Vertex_index{static_cast<unsigned>(vb->time_stamp() - 1)};
        auto [it, end] = CGAL::halfedges_around_source(index_va, mesh);
        // std::cerr << "  around mesh vertex " << index_va << ", search for vertex " << index_vb << '\n';
        // for(auto it2 = it; it2 != end; ++it2) {
        //   auto he = *it2;
        //   auto vd = target(he, mesh);
        //   std::cerr << "    " << vd << '\n';
        // }
        it = std::find_if(it, end, [&mesh, index_vb](auto he) { return target(he, mesh) == index_vb; });
        CGAL_assertion(it != end);
        auto he = CGAL::Euler::split_edge(*it, mesh);
        auto mesh_v = target(he, mesh);
        put(pmap, mesh_v, v->point());
        assert(mesh_v == Vertex_index{static_cast<unsigned>(time_stamp - 1)});
      }
      // for(auto e: edges(mesh)) {
      //   auto he = halfedge(e, mesh);
      //   auto vd1 = target(he, mesh);
      //   auto vd2 = source(he, mesh);
      //   if(!get(v_selected_map, vd1) || !get(v_selected_map, vd2)) continue;
      //   auto p1 = get(pmap, vd1);
      //   auto p2 = get(pmap, vd2);
      //   auto n = cdt.number_of_vertices();
      //   auto v1 = cdt.insert(p1);
      //   auto v2 = cdt.insert(p2);
      //   CGAL_assertion(n == cdt.number_of_vertices());
      //   auto steiner_vertices = cdt.sequence_of_Steiner_vertices(v1, v2);
      //   if(!steiner_vertices) continue;
      //   for(auto v: *steiner_vertices) {
      //     he = CGAL::Euler::split_edge(he, mesh);
      //     put(pmap, target(he, mesh), v->point());
      //   }
      // }
      std::ofstream out_mesh(options.dump_after_conforming_filename);
      out_mesh.precision(17);
      out_mesh << mesh;
      out_mesh.close();
    }

    if(!options.quiet) {
      std::cerr << "Number of vertices after conforming: " << cdt.number_of_vertices() << "\n\n";
    }
    CGAL_assertion(cdt.Delaunay::is_valid(true));
    CGAL_assertion(cdt.is_valid(true));
    CGAL_assertion(cdt.is_conforming());
    if(exit_code == EXIT_SUCCESS) {
      try {
        start_time = std::chrono::high_resolution_clock::now();
        cdt.restore_constrained_Delaunay();
        if(!options.quiet) {
          std::cout << "[timings] restored constrained Delaunay in " << std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::high_resolution_clock::now() - start_time).count() << " ms\n";
          std::cout << "Number of vertices after CDT: " << cdt.number_of_vertices() << "\n\n";
        }
      } catch(int error) {
        exit_code = error;
      }
    }
  } CDT_3_catch(CGAL::Failure_exception&) {
    finally();
    CDT_3_throw_exception_again;
  }
  finally();
  CGAL_assertion(cdt.is_conforming());
  CGAL_assertion(cdt.is_valid(true));


  return exit_code;
}
