#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel           Kernel;
typedef Kernel::Point_2                                             Point_2;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost      Cost;
typedef CGAL::Polyline_simplification_2::Stop_above_cost_threshold  Threshold;

int main(int, char**)
{
  std::vector<Point_2> points;

  points.push_back(Point_2(0, 0));
  points.push_back(Point_2(100, 0));
  points.push_back(Point_2(100, 100));
  points.push_back(Point_2(0, 100));
  points.push_back(Point_2(0, 50));
  points.push_back(Point_2(50, 50));
  points.push_back(Point_2(50, 51));
  points.push_back(Point_2(50, 49));
  points.push_back(Point_2(70, 49));

  std::vector<Point_2> output;
  CGAL::Polyline_simplification_2::simplify(points.begin(),
                                            points.end(),
                                            Cost(), Threshold (5.),
                                            std::back_inserter(output));
  return EXIT_SUCCESS;
}

