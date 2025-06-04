#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Random.h>
#include <CGAL/Memory_sizer.h>

#include <fstream>
#include <limits>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef int Color;

struct  Cell_data {
      unsigned char conflict_state;
};

typedef CGAL::Point_set_3<Point> Point_set;
typedef Point_set::Property_map<Color> Color_map;
typedef Point_set::Property_map<FT> FT_map;
typedef Point_set::Property_map<Cell_data> CD_map;
typedef Point_set::Property_map<bool> Bool_map;

typedef CGAL::Memory_sizer                                   Memory_sizer;


int main (int, char**)
{

  Memory_sizer memory_sizer;
  std::size_t base = memory_sizer.virtual_size();
  std::cout << "Memory usage: " << memory_sizer.virtual_size() << " bytes (virtual), "
            << memory_sizer.resident_size() << " bytes (resident)" << std::endl;

  Point_set point_set;

  Color black = 0;
  Color_map color2, color3, color4;
  CD_map cd_map;
  Bool_map bool_map;

  bool success;

  std::size_t  N = 1<<20; // 1 million points
  point_set.reserve (N); // For memory optimization
  for (std::size_t i = 0; i < N; ++ i)
  {
    Point_set::iterator it = point_set.insert (Point (double(i), double(i), double(i)));
  }

  double d;

  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per point " <<  std::endl;
  base = memory_sizer.virtual_size();

  std::tie (color2, success) = point_set.add_property_map<Color> ("color2", black);
  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();

  std::tie (color3, success) = point_set.add_property_map<Color> ("color3", black);
  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();


  std::tie (color4, success) = point_set.add_property_map<Color> ("color4", black);
  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();

  std::tie (cd_map, success) = point_set.add_property_map<Cell_data> ("cdm", Cell_data{0});

  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();

  std::tie (bool_map, success) = point_set.add_property_map<bool> ("bool", true);

  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();

  // Add more points to test memory usage
  for (std::size_t i = 0; i < 10; ++ i)
  {
    Point_set::iterator it = point_set.insert (Point (double(i), double(i), double(i)));
  }

  d= memory_sizer.virtual_size() - base;
  std::cout << "Memory usage: " << d / N << " bytes per property " <<  std::endl;
  base = memory_sizer.virtual_size();

  return 0;
}
