
#include <vector>
#include <CGAL/property_map.h>

int main()
{
  std::vector<int> v;
  CGAL::make_property_map(v);
  return 0;
}
