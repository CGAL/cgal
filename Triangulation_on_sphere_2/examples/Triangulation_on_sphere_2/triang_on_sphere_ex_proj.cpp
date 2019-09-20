#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/Projection_sphere_traits_3.h>

#include <boost/iterator/transform_iterator.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Projection_sphere_traits_3<K> Projection_traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Projection_traits> Projected_DToS2;

int main()
{
  std::vector<K::Point_3> points;
  points.push_back(K::Point_3( 3,  1, 1));
  points.push_back(K::Point_3(-8,  1, 1));
  points.push_back(K::Point_3( 1,  2, 1));
  points.push_back(K::Point_3( 1, -2, 1));

  Projection_traits traits(K::Point_3(1,1,1));
  Projected_DToS2 dtos(traits);

  Projection_traits::Construct_projected_point_3 cst =
    traits.construct_projected_point_3_object();

  for(const auto& pt : points)
  {
    std::cout << "Inserting: " << pt
              << " at distance: " << CGAL::squared_distance(pt, traits.center())
              << " from center" << std::endl;
    dtos.insert(cst(pt));

    std::cout << "dimension: " << dtos.dimension() << std::endl;
    std::cout << dtos.number_of_vertices() << " nv" << std::endl;
    std::cout << dtos.number_of_faces() << " nf" << std::endl;
    std::cout << dtos.number_of_ghost_faces() << " ng" << std::endl;
  }
}
