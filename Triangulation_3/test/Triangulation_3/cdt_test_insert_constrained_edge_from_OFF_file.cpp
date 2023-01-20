#define CGAL_DEBUG_CDT_3 1
#define CGAL_TRIANGULATION_CHECK_EXPENSIVE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Conforming_Delaunay_triangulation_3.h>
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

  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/fandisk.off");
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }


#ifndef CGAL_TEST_CDT_3_USE_CDT
  CGAL::Conforming_Delaunay_triangulation_3<Delaunay> cdt_edge;
#else
  CGAL::Constrained_Delaunay_triangulation_3<Delaunay> cdt_edge;
#endif
  auto point_map = get(CGAL::vertex_point, mesh);
  auto dt_vertex_handle_map =
      get(CGAL::dynamic_vertex_property_t<Delaunay::Vertex_handle>(), mesh);
  for(auto vertex_descriptor: vertices(mesh)) {
    auto vertex_handle = cdt_edge.insert(get(point_map, vertex_descriptor));
    put(dt_vertex_handle_map, vertex_descriptor, vertex_handle);
  }

  for(auto edge_descriptor: edges(mesh)) {
    auto s = source(edge_descriptor, mesh);
    auto t = target(edge_descriptor, mesh);
    cdt_edge.insert_constrained_edge(get(dt_vertex_handle_map, s),
                                     get(dt_vertex_handle_map, t));
  }
  {
    std::ofstream all_edges("all_segments.polylines.txt");
    all_edges.precision(17);
    cdt_edge.write_all_segments_file(all_edges);
  }
  assert(cdt_edge.is_conforming());

  return EXIT_SUCCESS;
}
