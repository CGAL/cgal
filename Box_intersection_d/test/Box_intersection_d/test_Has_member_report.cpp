#include <CGAL/assertions.h>
#include <CGAL/Box_intersection_d/segment_tree.h>
#include <cstddef>

struct S {
  void report(int);
  bool report();
};

struct With_report {
  bool report(int);
  bool report(int) const;
};

struct With_report_as_a_template_member_function {
  template <typename T> bool report(T);
};

int main() {
  using CGAL::Box_intersection_d::Has_member_report;
  CGAL_static_assertion(!Has_member_report<S>::value);
  CGAL_static_assertion(Has_member_report<With_report>::value);
  CGAL_static_assertion(Has_member_report<With_report_as_a_template_member_function>::value);
  return EXIT_SUCCESS;
}
