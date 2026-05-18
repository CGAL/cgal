// #define CGAL_MESH_2_VERBOSE 1
// #define CGAL_MESH_2_DEBUG_REFINEMENT_POINTS 1
// #define CGAL_MESHES_DEBUG_REFINEMENT_POINTS 1
// #define CGAL_MESH_2_DEBUG_BAD_EDGES 1
// #define CGAL_MESH_2_DEBUG_BAD_FACES 1
// #define CGAL_MESH_2_DEBUG_CLUSTERS 1
// #define CGAL_MESH_2_DEBUG_INSERTIONS 1
// #define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#include <CGAL/bisect_failures.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesher_no_edge_refinement_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/OBJ.h>
#include <CGAL/IO/write_VTU.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/Triangulation_simplex_base_with_time_stamp.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using FT = K::FT;
using Point_2 = K::Point_2;
using Vector_2 = K::Vector_2;

using Vb = CGAL::Triangulation_simplex_base_with_time_stamp<
    CGAL::Triangulation_vertex_base_with_info_2<std::size_t, K, CGAL::Triangulation_vertex_base_2<K>>>;
using Fb = CGAL::Delaunay_mesh_face_base_2<K>;
using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K, Tds, CGAL::Exact_predicates_tag>;

struct Obj_data {
  std::vector<Point_2> points;
  std::vector<std::vector<std::size_t> > polylines;
};

static auto save_obj(const Obj_data& data, const std::string& filename)
{
  struct Stream_plus_unused {
    std::ofstream out;
    std::vector<std::size_t> is_point_unused;
  } stream_plus_unused{ std::ofstream{filename}, {} };

  auto& out = stream_plus_unused.out;
  auto& is_point_unused = stream_plus_unused.is_point_unused;
  is_point_unused.resize(data.points.size(), 1);

  out.precision(17);
  if(!out)
    return stream_plus_unused;

  for(const std::vector<std::size_t>& pl : data.polylines) {
    if(pl.size() < 2)
      continue;
    for(std::size_t id : pl) {
      is_point_unused[id] = 0;
    }
  }

  for(const Point_2& p : data.points) {
    if(is_point_unused[&p - &data.points[0]] == 0) {
      out << "v " << p.x() << " " << p.y() << " 0\n";
    }
  }

  std::partial_sum(is_point_unused.begin(), is_point_unused.end(), is_point_unused.begin());
  for(const std::vector<std::size_t>& pl : data.polylines) {
    if(pl.size() < 2)
      continue;
    out << "l";
    for(std::size_t id : pl) {
      out << " " << (id + 1 - is_point_unused[id]);
    }
    out << "\n";
  }
  return stream_plus_unused;
}

static int run_mesh(const Obj_data& data)
{
  CDT cdt;

  double squared_distance = (std::numeric_limits<double>::max)();
  for(const std::vector<std::size_t>& id_pl : data.polylines) {
    for(std::size_t pid = 1; pid < id_pl.size(); ++pid) {
      assert(data.points[id_pl[pid - 1]] != data.points[id_pl[pid]]);
      squared_distance =
          (std::min)(squared_distance, CGAL::squared_distance(data.points[id_pl[pid - 1]], data.points[id_pl[pid]]));
      auto va = cdt.insert(data.points[id_pl[pid - 1]]);
      va->info() = id_pl[pid - 1];
      auto vb = cdt.insert(data.points[id_pl[pid]]);
      vb->info() = id_pl[pid];
      // std::cerr << "Inserting constraint between " << CGAL::IO::oformat(va, CGAL::With_point_tag{})
      //           << " and " << CGAL::IO::oformat(vb, CGAL::With_point_tag{}) << "\n";
      cdt.insert_constraint(va, vb);
    }
  }

  auto [file_stream, points_ids_offsets] = save_obj(data, "initial_cdt.obj");

  if(cdt.dimension() == 2) {
    CDT::Face_circulator fc = cdt.incident_faces(cdt.infinite_vertex()), done = fc;
    do{
      auto edge_index = fc->index(cdt.infinite_vertex());
      if(!fc->is_constrained(edge_index)) {
        auto va = fc->vertex(cdt.cw(edge_index));
        auto vb = fc->vertex(cdt.ccw(edge_index));
      //   std::cerr << "Add constrained edge on the convex hull: ";
      //   std::cerr << CGAL::IO::oformat(va, CGAL::With_point_tag{}) << "  "
      //             << CGAL::IO::oformat(vb, CGAL::With_point_tag{}) << "\n";
        cdt.insert_constraint(va, vb);
        file_stream << "l " << (va->info() + 1 - points_ids_offsets[va->info()]) << " "
                    << (vb->info() + 1 - points_ids_offsets[vb->info()]) << "\n";
      }
      ++fc;
    } while (fc != done);
  }

  file_stream.close();

  std::cout << "Before conforming Gabriel: "
            << cdt.number_of_vertices() << " vertices.\n"
            << "  Smallest squared distance between constraint endpoints: "
            << squared_distance << "\n";

  // CGAL::refine_Delaunay_mesh_2(cdt, CGAL::parameters::criteria(CGAL::Delaunay_mesh_criteria_2<CDT>{}));

  CGAL::Triangulation_conformer_2<CDT> conform(cdt);
  conform.init_Gabriel();
  // auto i = 0u;
  while(!conform.is_conforming_done()) {
    conform.try_one_step_conforming_Gabriel();
    // std::string filename = "debug-conforming-gabriel-step-" + std::to_string(i) + ".vtu";
    // CGAL::IO::write_VTU(filename, cdt, CGAL::IO::ASCII);
    // ++i;
  }
  CGAL_assertion(cdt.is_valid(true));
  std::cout << "done" << std::endl;
  CGAL::IO::write_VTU("conformed_cdt.vtu", cdt, CGAL::IO::ASCII);
  return EXIT_SUCCESS;
}

static std::string join_path(const std::string& dir, const std::string& name)
{
  if(dir.empty())
    return name;
  if(dir.back() == '/')
    return dir + name;
  return dir + "/" + name;
}

int main(int argc, char*argv[] )
{
  std::cerr.precision(17);
  std::cout.precision(17);
  std::clog.precision(17);
  std::ifstream in(argc > 1 ? argv[1] : "mini.obj");
  if(!in)
    return EXIT_FAILURE;

  Obj_data data;
  std::vector<std::vector<std::size_t> > unused_id_polygons;
  bool success = CGAL::IO::internal::read_OBJ(in, data.points, data.polylines, unused_id_polygons);
  if(!success)
    return EXIT_FAILURE;

  auto get_size_fn = [](const Obj_data& d) -> std::size_t {
    return d.polylines.size();
  };

  auto simplify_fn = [](Obj_data& d, std::size_t start, std::size_t end) -> bool {
    if(start >= end || start >= d.polylines.size())
      return false;
    end = (std::min)(end, d.polylines.size());
    d.polylines.erase(d.polylines.begin() + start, d.polylines.begin() + end);
    return true;
  };

  auto run_fn = [](const Obj_data& d) -> int {
    return run_mesh(d);
  };

  std::string output_dir = "./";

  auto save_fn = [&](const Obj_data& d, const std::string& prefix) {
    const std::string filename = join_path(output_dir, "bisect_" + prefix + ".obj");
    save_obj(d, filename);
  };

  if(argc > 2) {
    output_dir = argv[2];
    return CGAL::bisect_failures(data, get_size_fn, simplify_fn, run_fn, save_fn);
  } else {
    return run_mesh(data);
  }
}
