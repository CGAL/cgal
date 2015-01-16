
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/AABB_halfedge_graph_segment_primitive.h>

#include <CGAL/Polygon_mesh_slicer_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;

int main(int, char** argv)
{
  std::ifstream input("data/U.off");
  Mesh m;

  if (!input || !(input >> m)){
    std::cerr << "Error: can not read file.";
    return 1;
  }

  AABB_tree tree(edges(m).first, edges(m).second, m);

  CGAL::Polygon_mesh_slicer_3<Mesh, K> slicer(m, tree);
  
  return 0;
}
