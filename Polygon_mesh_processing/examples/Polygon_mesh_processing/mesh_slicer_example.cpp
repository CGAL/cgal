#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_slicer.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/AABB_halfedge_graph_segment_primitive.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <iostream>
#include <list>
#include <string>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

typedef std::vector<K::Point_3>                               Polyline_type;
typedef std::list<Polyline_type>                              Polylines;

typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh>     HGSP;
typedef CGAL::AABB_traits<K, HGSP>                            AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                          AABB_tree;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh) || CGAL::is_empty(mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  // Slicer constructor from the mesh
  CGAL::Polygon_mesh_slicer<Mesh, K> slicer(mesh);

  Polylines polylines;
  slicer(K::Plane_3(0, 0, 1, -0.4), std::back_inserter(polylines));
  std::cout << "At z = 0.4, the slicer intersects "
            << polylines.size() << " polylines" << std::endl;
  polylines.clear();

  slicer(K::Plane_3(0, 0, 1, 0.2), std::back_inserter(polylines));
  std::cout << "At z = -0.2, the slicer intersects "
            << polylines.size() << " polylines" << std::endl;
  polylines.clear();

  // Use the Slicer constructor from a pre-built AABB_tree
  AABB_tree tree(edges(mesh).first, edges(mesh).second, mesh);
  CGAL::Polygon_mesh_slicer<Mesh, K> slicer_aabb(mesh, tree);
  slicer_aabb(K::Plane_3(0, 0, 1, -0.4), std::back_inserter(polylines));
  std::cout << "At z = 0.4, the slicer intersects "
            << polylines.size() << " polylines" << std::endl;
  polylines.clear();

  return 0;
}
