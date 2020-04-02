#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Regular_triangulation.h>
#include <CGAL/assertions.h>

#include <iostream>
#include <iterator>
#include <vector>

const int D = 5;    // Dimension
const int N = 100;  // Number of points
typedef CGAL::Epick_d< CGAL::Dimension_tag<D> >         K;
typedef CGAL::Regular_triangulation<K>                  T;
typedef K::Point_d                                      Bare_point;
typedef K::Weighted_point_d                             Weighted_point;

int main()
{
  // Instantiate a random point generator
  CGAL::Random rng(0);
  typedef CGAL::Random_points_in_cube_d<Bare_point> Random_points_iterator;
  Random_points_iterator rand_it(D, 1.0, rng);

  // Generate N random points
  std::vector<Weighted_point> points;
  for (int i = 0; i < N; ++i)
    points.push_back(Weighted_point(*rand_it++, rng.get_double(0., 10.)));

  T t(D);
  CGAL_assertion(t.empty());

  // Insert the points in the triangulation
  t.insert(points.begin(), points.end());
  CGAL_assertion( t.is_valid() );
  std::cout << "Regular triangulation successfully computed: "
    << t.number_of_vertices() << " vertices, "
    << t.number_of_finite_full_cells() << " finite cells."
    << std::endl;

  return 0;
}
