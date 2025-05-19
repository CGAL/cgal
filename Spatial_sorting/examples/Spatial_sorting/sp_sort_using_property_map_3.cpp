#include <CGAL/Simple_cartesian.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <vector>
#include <CGAL/boost/iterator/counting_iterator.hpp>

typedef CGAL::Simple_cartesian<double>                  Kernel;
typedef Kernel::Point_3                                 Point_3;
typedef CGAL::Spatial_sort_traits_adapter_3<Kernel,
          CGAL::Pointer_property_map<Point_3>::type > Search_traits_3;

int main()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(1,3,11));
  points.push_back(Point_3(14,34,46));
  points.push_back(Point_3(414,34,4));
  points.push_back(Point_3(4,2,56));
  points.push_back(Point_3(744,4154,43));
  points.push_back(Point_3(74,44,1));

  std::vector<std::size_t> indices;
  indices.reserve(points.size());

  std::copy(boost::counting_iterator<std::size_t>(0),
            boost::counting_iterator<std::size_t>(points.size()),
            std::back_inserter(indices));

  std::cout << "Order using default policy (median)\n";
  CGAL::spatial_sort( indices.begin(),
                      indices.end(),
                      Search_traits_3(CGAL::make_property_map(points)) );

  for (std::size_t i : indices)
    std::cout << points[i] << "\n";

  std::cout << "Order using middle policy\n";
  CGAL::spatial_sort( indices.begin(),
                      indices.end(),
                      Search_traits_3(CGAL::make_property_map(points)),
                      CGAL::Hilbert_sort_middle_policy());

  for (std::size_t i : indices)
    std::cout << points[i] << "\n";

  std::cout << "done" << std::endl;

  return 0;
}
