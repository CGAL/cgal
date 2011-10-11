#include <CGAL/Simple_cartesian.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <vector>
#include <boost/iterator/counting_iterator.hpp>

typedef CGAL::Simple_cartesian<double>                  Kernel;
typedef Kernel::Point_3                                 Point_3;
//using a pointer as a special property map type
typedef 
  CGAL::Spatial_sort_traits_adapter_3<Kernel,Point_3*>  Search_traits_3;

int main()
{
  std::vector<Point_3> points;
  points.push_back(Point_3(1,3,11));
  points.push_back(Point_3(14,34,46));
  points.push_back(Point_3(414,34,4));
  points.push_back(Point_3(4,2,56));
  points.push_back(Point_3(744,4154,43));
  points.push_back(Point_3(74,44,1));
  
  std::vector<std::ptrdiff_t> indices;
  indices.reserve(points.size());
  
  std::copy(boost::counting_iterator<std::ptrdiff_t>(0),
            boost::counting_iterator<std::ptrdiff_t>(points.size()),
            std::back_inserter(indices));
  
  CGAL::spatial_sort( indices.begin(),indices.end(),Search_traits_3(&(points[0])) );

  for (std::vector<std::ptrdiff_t>::iterator it=indices.begin();it!=indices.end();++it)
    std::cout << points[*it] << "\n";

  std::cout << "done" << std::endl;
  
  return 0;
}
