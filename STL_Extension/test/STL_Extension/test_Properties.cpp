
#include <CGAL/Properties.h>

using namespace CGAL::Properties;

int main() {

  Property_container properties;

  auto &a = properties.get<int, decltype("A"_param)>();
  auto &b = properties.get<int, decltype("B"_param)>();
  auto &c = properties.get<float, decltype("C"_param)>();

  assert((std::is_same<decltype("A"_param), decltype("A"_param)>::value));

  assert(a != b);

  auto &a2 = properties.get<int, decltype("A"_param)>();
  assert(a == a2);

  properties.threadsafe_get<float, decltype("C"_param)>();

  return 0;
}