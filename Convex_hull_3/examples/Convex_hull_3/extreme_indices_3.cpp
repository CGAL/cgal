#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Extreme_points_traits_adapter_3.h>
#include <CGAL/IO/read_off_points.h>

#include <boost/iterator/counting_iterator.hpp>

#include <vector>
#include <fstream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel      K;
typedef K::Point_3                                               Point_3;

int main(int argc, char* argv[])
{
  std::ifstream in( (argc>1)? argv[1] : "data/star.off");

  if(!in)
  {
    std::cerr<< "Cannot open input file." <<std::endl;
    return 1;
  }

  std::vector<Point_3> points;
  CGAL::read_off_points(in, std::back_inserter(points));

  //This will contain the extreme vertices
  std::vector<std::size_t> extreme_point_indices;

  //call the function with the traits adapter for vertices
  CGAL::extreme_points_3(CGAL::make_range(boost::counting_iterator<std::size_t>(0),
                                          boost::counting_iterator<std::size_t>(points.size())),
                         std::back_inserter(extreme_point_indices),
                         CGAL::make_extreme_points_traits_adapter(CGAL::make_property_map(points)));
  //print the number of extreme vertices
  std::cout << "Indices of points on the convex hull are:\n";
  std::copy(extreme_point_indices.begin(), extreme_point_indices.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << "\n";

  return 0;
}
