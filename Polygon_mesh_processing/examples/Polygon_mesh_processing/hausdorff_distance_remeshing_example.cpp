#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#define TAG CGAL::Parallel_if_available_tag

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3                     Point;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
  Mesh tm1, tm2;
  CGAL::make_tetrahedron(Point(.0,.0,.0),
                         Point(2,.0,.0),
                         Point(1,1,1),
                         Point(1,.0,2),
                         tm1);
  tm2=tm1;
  CGAL::Polygon_mesh_processing::isotropic_remeshing(tm2.faces(),.05, tm2);

  std::cout << "Approximated Hausdorff distance: "
            << CGAL::Polygon_mesh_processing::approximate_Hausdorff_distance
                  <TAG>(tm1, tm2, PMP::parameters::number_of_points_per_area_unit(4000))
            << std::endl;
}
