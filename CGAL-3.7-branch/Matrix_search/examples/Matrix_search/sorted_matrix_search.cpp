#include <CGAL/Random.h>
#include <CGAL/Cartesian_matrix.h>
#include <CGAL/sorted_matrix_search.h>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>

typedef int                                     Value;
typedef std::vector<Value>                      Vector;
typedef Vector::iterator                        Value_iterator;
typedef std::vector<Vector>                     Vector_cont;
typedef CGAL::Cartesian_matrix<std::plus<int>,
                               Value_iterator,
                               Value_iterator>  Matrix;

int main()
{
  // set of vectors the matrices are build from:
  Vector_cont vectors;

  // generate a random vector and sort it:
  Vector a;
  const int n = 5;
  for (int i = 0; i < n; ++i)
    a.push_back(CGAL::default_random(100));
  std::sort(a.begin(), a.end());
  std::cout << "a = ( ";
  std::copy(a.begin(), a.end(), std::ostream_iterator<int>(std::cout," "));
  std::cout << ")\n";

  // build a Cartesian matrix from a:
  Matrix M(a.begin(), a.end(), a.begin(), a.end());

  // search for an upper bound for max(a):
  Value bound = a[n-1];
  Value upper_bound =
  CGAL::sorted_matrix_search(
    &M, &M + 1,
    CGAL::sorted_matrix_search_traits_adaptor(
      std::bind2nd(std::greater_equal<Value>(), bound), M));
  std::cout << "Upper bound for " << bound << " is "
            << upper_bound << "." << std::endl;

  return 0;
}
