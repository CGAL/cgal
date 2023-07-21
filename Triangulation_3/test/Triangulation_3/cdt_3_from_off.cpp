//#define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#define NO_TRY_CATCH 1
#define CGAL_DEBUG_CDT_3 1
#define CGAL_TRIANGULATION_CHECK_EXPENSIVE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/Base_with_time_stamp.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <CGAL/Polygon_mesh_processing/bbox.h>

#include <vector>
#include <cassert>
#include <fstream>
#include <string>

#if NO_TRY_CATCH
// Iff -fno-exceptions, transform error handling code to work without it.
# define CDT_3_try      if (true)
# define CDT_3_catch(X) if (false)
# define CDT_3_throw_exception_again
#else
// Else proceed normally.
# define CDT_3_try      try
# define CDT_3_catch(X) catch(X)
# define CDT_3_throw_exception_again throw
#endif

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb = CGAL::Base_with_time_stamp<CGAL::Constrained_Delaunay_triangulation_vertex_base_3<K>>;
using Cb = CGAL::Constrained_Delaunay_triangulation_cell_base_3<K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using Delaunay = CGAL::Delaunay_triangulation_3<K, Tds>;
using Point = Delaunay::Point;
using Point_3 = K::Point_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor =boost::graph_traits<Mesh>::face_descriptor;

int main(int argc, char* argv[])
{
  std::cerr.precision(17);
  std::cout.precision(17);

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/mpi.off");
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  const std::string output_filename = (argc > 2) ? argv[2] : "dump.off";

  CGAL::Constrained_Delaunay_triangulation_3<Delaunay> cdt;

  const auto bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  double d_x = bbox.xmax() - bbox.xmin();
  double d_y = bbox.ymax() - bbox.ymin();
  double d_z = bbox.zmax() - bbox.zmin();

  const double max_d = (std::max)(d_x, (std::max)(d_y, d_z));
  if(d_x == 0) d_x = max_d;
  if(d_y == 0) d_y = max_d;
  if(d_z == 0) d_z = max_d;

  // cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
  // cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
  // cdt.insert(Point(bbox.xmin() - d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
  // cdt.insert(Point(bbox.xmin() - d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));
  // cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmin() - d_z));
  // cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmin() - d_z));
  // cdt.insert(Point(bbox.xmax() + d_x, bbox.ymin() - d_y, bbox.zmax() + d_z));
  // cdt.insert(Point(bbox.xmax() + d_x, bbox.ymax() + d_y, bbox.zmax() + d_z));

  int exit_code = EXIT_SUCCESS;

  auto finally = [&cdt,output_filename]() {
    {
      std::ofstream dump("dump.binary.cgal");
      CGAL::IO::save_binary_file(dump, cdt);
    }
    {
      std::ofstream dump(output_filename);
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

  auto pmap = get(CGAL::vertex_point, mesh);
  for(auto v: vertices(mesh)) {
    cdt.insert(get(pmap, v));
  }
  int poly_id = 0;
  CDT_3_try {
    for(auto face_descriptor : faces(mesh)) {
      std::vector<Point_3> polygon;
      const auto he = halfedge(face_descriptor, mesh);
      for(auto vertex_it : CGAL::vertices_around_face(he, mesh)) {
        polygon.push_back(get(pmap, vertex_it));
      }
      std::cerr << "NEW POLYGON #" << poly_id << '\n';
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
        auto id = cdt.insert_constrained_polygon(polygon, true);
        assert(id == poly_id);
        ++poly_id;
      } catch(int error) {
        exit_code = error;
      }
      // std::ofstream dump("dump.binary.cgal");
      // CGAL::Mesh_3::save_binary_file(dump, cdt);
    }
    cdt.restore_Delaunay();
    for(auto e: edges(mesh)) {
      auto he = halfedge(e, mesh);
      auto p1 = get(pmap, target(he, mesh));
      auto p2 = get(pmap, source(he, mesh));
      auto n = cdt.number_of_vertices();
      auto v1 = cdt.insert(p1);
      auto v2 = cdt.insert(p2);
      CGAL_assertion(n == cdt.number_of_vertices());  
      auto steiner_vertices = cdt.sequence_of_Steiner_vertices(v1, v2);
      for(auto v: steiner_vertices) {
        he = CGAL::Euler::split_edge(he, mesh);
        put(pmap, target(he, mesh), v->point());
      }
    }
    std::ofstream out_mesh("out-conforming.off");
    out_mesh.precision(17);
    out_mesh << mesh;
    out_mesh.close();
    {
      std::ofstream all_edges("dump_all_segments.polylines.txt");
      all_edges.precision(17);
      cdt.write_all_segments_file(all_edges);
    }

    std::cerr << "Number of vertices after conforming: " << cdt.number_of_vertices() << '\n';
    assert(cdt.is_conforming());
    if(exit_code == EXIT_SUCCESS) {
      try {
        cdt.restore_constrained_Delaunay();
        std::cerr << "Number of vertices after CDT: " << cdt.number_of_vertices() << '\n';
      } catch(int error) {
        exit_code = error;
      }
    }
  } CDT_3_catch(CGAL::Failure_exception&) {
    finally();
    CDT_3_throw_exception_again;
  }
  finally();
  assert(cdt.is_conforming());

  return exit_code;
}
