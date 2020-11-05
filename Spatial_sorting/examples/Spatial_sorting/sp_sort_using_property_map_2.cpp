#include <CGAL/Simple_cartesian.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_2.h>
#include <CGAL/property_map.h>
#include <vector>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef std::pair<Point_2,int>              Point_with_info;
typedef std::vector< Point_with_info >      Data_vector;

typedef CGAL::Spatial_sort_traits_adapter_2<Kernel,
  CGAL::First_of_pair_property_map<Point_with_info>
> Search_traits_2;

int main()
{
  Data_vector points;
  points.push_back(std::make_pair(Point_2(14,12) , 3));
  points.push_back(std::make_pair(Point_2(1,2)   , 0));
  points.push_back(std::make_pair(Point_2(414,2) , 5));
  points.push_back(std::make_pair(Point_2(4,21)  , 1));
  points.push_back(std::make_pair(Point_2(7,74)  , 2));
  points.push_back(std::make_pair(Point_2(74,4)  , 4));

  Search_traits_2 traits;
  CGAL::spatial_sort(points.begin(), points.end(), traits);
  for (Data_vector::iterator it=points.begin();it!=points.end();++it)
    std::cout << it->second << " ";
  std::cout << "\n";

  std::cout << "done" << std::endl;

  return 0;
}
