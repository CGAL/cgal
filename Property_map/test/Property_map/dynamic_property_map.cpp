#include <CGAL/Dynamic_property_map.h>

int main()
{
  CGAL::internal::Dynamic_property_map<int,int> dpm(2);
  assert(dpm.default_value() == 2);

  assert(get(dpm, 0) == 2);
  put(dpm, 0, 1);
  assert(get(dpm, 0) == 1);

  return 0;
}
