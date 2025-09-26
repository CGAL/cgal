#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel   Epick;

int main()
{
CGAL::Surface_mesh<CGAL::Epick::Point_3> mesh;
  auto f0 = mesh.add_face();
  auto f1 = mesh.add_face();  // If this line is commented out, it will not crash.
  mesh.remove_face(f0);

  auto v0 = mesh.add_vertex();
  auto v1 = mesh.add_vertex();  // If this line is commented out, it will not crash.
  mesh.remove_vertex(v0);

  auto e0 = mesh.add_edge();
  auto e1 = mesh.add_edge();  // If this line is commented out, it will not crash.
  mesh.remove_edge(edge(e0,mesh));


  mesh.collect_garbage();

  return 0;
}