#define CGAL_DEBUG_CDT_3 1
#define CGAL_TRIANGULATION_CHECK_EXPENSIVE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <vector>
#include <cassert>
#include <fstream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     K;
typedef CGAL::Triangulation_data_structure_3<
  CGAL::Constrained_Delaunay_triangulation_vertex_base_3<K>,
  CGAL::Constrained_Delaunay_triangulation_cell_base_3<K> >     Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds>                  Delaunay;
typedef Delaunay::Point                                         Point;
using Point_3 = K::Point_3;

typedef CGAL::Surface_mesh<Point>                      Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor          face_descriptor;

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

  auto pmap = get(CGAL::vertex_point, mesh);
  for(auto face_descriptor: faces(mesh)) {
    std::vector<Point_3> polygon;
    const auto he = halfedge(face_descriptor, mesh);
    for(auto vertex_it: CGAL::vertices_around_face(he, mesh)) {
      polygon.push_back(get(pmap, vertex_it));
    }
    std::cerr << "NEW POLYGON\n";
    try {
      cdt.insert_constrained_polygon(polygon);
    } catch (int error) {
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
  assert(cdt.is_conforming());

  return exit_code;
}
