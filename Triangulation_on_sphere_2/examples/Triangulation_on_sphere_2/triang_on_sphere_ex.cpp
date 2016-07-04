#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_sphere_2.h>
#include <CGAL/Delaunay_triangulation_sphere_traits_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Delaunay_triangulation_sphere_traits_2<K> Traits;
typedef CGAL::Delaunay_triangulation_sphere_2<Traits> DToS2;

int main()
{
  std::vector<K::Point_3> points;
  points.push_back( K::Point_3(2,1,1) );
  points.push_back( K::Point_3(-2,1,1) );
  points.push_back( K::Point_3(1,2,1) );
  points.push_back( K::Point_3(1,-2,1) );


  Traits traits(K::Point_3(1,1,1));
  DToS2 dtos(traits);
  dtos.insert(points.begin(),points.end());
}


