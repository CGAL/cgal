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

  CGAL::Constrained_Delaunay_triangulation_3<Delaunay> cdt;

  int exit_code = EXIT_SUCCESS;

  auto finally = [&cdt]() {
    {
      std::ofstream dump("dump.binary.cgal");
      CGAL::Mesh_3::save_binary_file(dump, cdt);
    }
    {
      std::ofstream missing_faces("missing_faces.polylines.txt");
      missing_faces.precision(17);
      cdt.write_missing_subfaces_file(missing_faces);
    }
    {
      std::ofstream missing_edges("missing_segments.polylines.txt");
      missing_edges.precision(17);
      if(cdt.write_missing_segments_file(missing_edges)) {
        std::cerr << "ERROR: Missing segments!\n";
      }
    }
  };

  auto pmap = get(CGAL::vertex_point, mesh);
  int poly_id = 0;
  CDT_3_try {
    for(auto face_descriptor : faces(mesh)) {
      std::vector<Point_3> polygon;
      const auto he = halfedge(face_descriptor, mesh);
      for(auto vertex_it : CGAL::vertices_around_face(he, mesh)) {
        polygon.push_back(get(pmap, vertex_it));
      }
      std::cerr << "NEW POLYGON #" << poly_id << '\n';
      try {
        auto id = cdt.insert_constrained_polygon(polygon);
        assert(id == poly_id);
        ++poly_id;
      } catch(int error) {
        exit_code = error;
      }
      // std::ofstream dump("dump.binary.cgal");
      // CGAL::Mesh_3::save_binary_file(dump, cdt);
    }
    assert(cdt.is_conforming());
    if(exit_code == EXIT_SUCCESS) {
      try {
        cdt.restore_constrained_Delaunay();
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
