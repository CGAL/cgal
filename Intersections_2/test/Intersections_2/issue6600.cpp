#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>

int main()
{
  typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
  typedef Kernel::Point_3 Point;
  typedef CGAL::Surface_mesh<Point> Surface_mesh;
  Surface_mesh mesh;
  CGAL::Polygon_mesh_processing::experimental::autorefine_and_remove_self_intersections(mesh);

  return 0;
}
