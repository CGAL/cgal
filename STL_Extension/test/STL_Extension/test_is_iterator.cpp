#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <CGAL/is_iterator.h>
#include <vector>

int main() {
  typedef std::vector<int>::const_iterator vector_it;
  typedef int* int_p;
  using CGAL::is_iterator;

  CGAL_static_assertion(is_iterator<vector_it>::value);
  CGAL_static_assertion(!is_iterator<void>::value);
  CGAL_static_assertion(!is_iterator<int>::value);
  CGAL_static_assertion(is_iterator<int_p>::value);
}
