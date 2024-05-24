#include <iostream>
#include <map>

#include <CGAL/config.h>
#if __has_include(<source_location>)
#include <CGAL/NT_wrapper.h>

struct Map {
  std::map<std::string, std::size_t> m;
  void increment(const std::string& key) {
    m[key]++;
  }
  ~Map() {
    for(const auto& [key, value] : m)
      std::clog << value << " " << key << std::endl;
  }
} map;

void incr_counter(const std::source_location& l) {
  map.increment(l.function_name());
}

using K = CGAL::Wrapped_epick;

using EK = CGAL::Wrapped_epeck;

int main() {
  CGAL::NT_wrapper<double>::f = incr_counter;
  CGAL::NT_wrapper<CGAL::Epeck_ft>::f = incr_counter;
  K ::Point_2 p1(0, 2), p2(2, 4), p3(4, 6);
  EK::Point_2 p4(0, 2), p5(2, 4), p6(4, 6);
  exact(p4);
  exact(p5);
  exact(p6);
  [[maybe_unused]] bool result_2d_k = CGAL::orientation(p1, p2, p3) == CGAL::COLLINEAR;
  [[maybe_unused]] bool result_2d_e = CGAL::orientation(p4, p5, p6) == CGAL::COLLINEAR;
  assert (result_2d_k && result_2d_e);

  for(auto [key, _] : map.m) {
    if(key.find("const") != std::string::npos) {
      std::cerr << "ERROR: copy constructor called!\n";
      return 1;
    }
  }
  return 0;
}

#else
int main() {
  std::cerr << "C++20 <source_location> not available\n";
  return 0;
}
#endif
