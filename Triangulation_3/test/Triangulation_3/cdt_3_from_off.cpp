#define CGAL_DEBUG_CDT_3 1
#define CGAL_TRIANGULATION_CHECK_EXPENSIVE 1
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Random.h>
#include <CGAL/Constrained_Delaunay_triangulation_3.h>
#include <CGAL/draw_triangulation_3.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/IO/File_binary_mesh_3.h>

#include <vector>
#include <cassert>
#include <fstream>

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

  const char* filename = (argc > 1) ? argv[1] : "data/fandisk.off";
  std::ifstream input(filename);
  Mesh mesh;
  if (!input || !(input >> mesh))
  {
    std::cerr << "Not a valid input file." << std::endl;
    return 1;
  }

  CGAL::Constrained_Delaunay_triangulation_3<Delaunay> cdt;

  auto pmap = get(CGAL::vertex_point, mesh);
  for(auto face_descriptor: faces(mesh)) {
    std::vector<Point_3> polygon;
    const auto he = halfedge(face_descriptor, mesh);
    for(auto vertex_it: CGAL::vertices_around_face(he, mesh)) {
      polygon.push_back(get(pmap, vertex_it));
    }
    std::cerr << "NEW POLYGON\n";
    cdt.insert_constrained_polygon(polygon);
    std::ofstream dump("dump.binary.cgal");
    CGAL::Mesh_3::save_binary_file(dump, cdt);
  }

  // auto edge_is_feature_map = get(CGAL::edge_is_feature, mesh);
  // auto face_patch_id_map = get(CGAL::face_patch_id_t<int>(), mesh);

  // namespace PMP = CGAL::Polygon_mesh_processing;
  // PMP::sharp_edges_segmentation(mesh, 80, edge_is_feature_map,
  //                               face_patch_id_map);

  // auto point_map = get(CGAL::vertex_point, mesh);
  // auto dt_vertex_handle_map =
  //     get(CGAL::dynamic_vertex_property_t<Delaunay::Vertex_handle>(), mesh);
  // for(auto vertex_descriptor: vertices(mesh)) {
  //   auto vertex_handle = cdt.insert(get(point_map, vertex_descriptor));
  //   put(dt_vertex_handle_map, vertex_descriptor, vertex_handle);
  // }

  // for(auto edge_descriptor: edges(mesh)) {
  //   if(get(edge_is_feature_map, edge_descriptor)) {
  //     auto s = source(edge_descriptor, mesh);
  //     auto t = target(edge_descriptor, mesh);
  //     cdt.insert_constrained_edge(get(dt_vertex_handle_map, s),
  //                                get(dt_vertex_handle_map, t));
  //   }
  // }
  CGAL::draw(static_cast<Delaunay&>(cdt), "CDT_3", true);
  assert(cdt.is_conforming());
}
