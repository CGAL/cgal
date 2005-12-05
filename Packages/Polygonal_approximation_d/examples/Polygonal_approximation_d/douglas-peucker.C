//! \file examples/Polygonal_approximation_d/douglas-peucker.C
#include <CGAL/Simple_cartesian.h>
#include <CGAL/simplify_polyline.h>
#include <list>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_2 Point_2;
typedef CGAL::Squared_euclidean_error<K,CGAL::Max_tag,CGAL::Segment_tag,2> Error;

int main()
{
  std::list<Point_2> polyline, result;
  double eps;
  std::size_t n;

  std::cin >> n;
  std::copy(std::istream_iterator<Point_2>(std::cin), 
	    std::istream_iterator<Point_2>(), 
	    std::back_inserter(polyline));

  CGAL::simplify_polyline_bound_number_of_points<CGAL::Recursive_split, CGAL::Global, Error>
    (polyline.begin(), polyline.end(),
     n/2,  eps,
     std::back_inserter(result));
  
  std::cerr << "Maximal error: " << eps << std::endl;
  std::cout << n/2 << std::endl;
  std::copy(result.begin(), result.end(), std::ostream_iterator<Point_2>(std::cout, "\n"));
  
  return 0;
}
