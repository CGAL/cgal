#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/regular_neighbor_coordinates_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Regular_triangulation_euclidean_traits_2<K> Gt;
typedef CGAL::Regular_triangulation_2<Gt>              Regular_triangulation;
typedef Regular_triangulation::Weighted_point Weighted_point;
typedef std::vector< std::pair< Weighted_point, K::FT  > >
                                                       Point_coordinate_vector;

int main()
{
  Regular_triangulation rt;

  for (int y=0 ; y<3 ; y++)
    for (int x=0 ; x<3 ; x++)
      rt.insert(Weighted_point(K::Point_2(x,y), 0));

  //coordinate computation
  Weighted_point wp(K::Point_2(1.2, 0.7),2);
  Point_coordinate_vector  coords;
  CGAL::Triple<
    std::back_insert_iterator<Point_coordinate_vector>,
    K::FT, bool> result =
    CGAL::regular_neighbor_coordinates_2(rt, wp,
					 std::back_inserter(coords));
  if(!result.third){
    std::cout << "The coordinate computation was not successful."
	      << std::endl;
    std::cout << "The point (" <<wp.point() << ") lies outside the convex hull."
	      << std::endl;
  }
  K::FT  norm = result.second;
  std::cout << "Coordinate computation successful." << std::endl;
  std::cout << "Normalization factor: " <<norm << std::endl;

  std::cout << "done" << std::endl;
  return 0;
}
