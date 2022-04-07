#include <CGAL/assertions.h>
#include <CGAL/is_iterator.h>
#include <CGAL/use.h>
#include <vector>
#include <list>

struct A { };

int main() {
  typedef std::vector<int>::const_iterator vector_it;
  typedef std::list<int>::const_iterator list_it;
  typedef int* int_p;
  using CGAL::is_iterator;
  using CGAL::is_iterator_type;
  using CGAL::is_iterator_to;

  CGAL_USE_TYPE(vector_it);
  CGAL_USE_TYPE(list_it);
  CGAL_USE_TYPE(int_p);
  CGAL_static_assertion(is_iterator<vector_it>::value);
  CGAL_static_assertion(is_iterator<list_it>::value);
  CGAL_static_assertion(!is_iterator<void>::value);
  CGAL_static_assertion(!is_iterator<int>::value);
  CGAL_static_assertion(is_iterator<int_p>::value);

  CGAL_static_assertion((is_iterator_type<vector_it,std::bidirectional_iterator_tag>::value));
  CGAL_static_assertion((!is_iterator_type<list_it,std::random_access_iterator_tag>::value));
  CGAL_static_assertion((!is_iterator_type<short,std::output_iterator_tag>::value));

  CGAL_static_assertion((is_iterator_to<int_p,double>::value));
  CGAL_static_assertion((!is_iterator_to<A,int>::value));
  CGAL_static_assertion((!is_iterator_to<A*,int>::value));
}
