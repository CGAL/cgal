#include <CGAL/Random.h>
#include <CGAL/function_objects.h>
#include <CGAL/Cartesian_matrix.h>
#include <CGAL/sorted_matrix_search.h>
#include <vector>

using namespace std;
using namespace CGAL;

typedef int                              Value;
typedef vector< Value >                  Vector;
typedef Vector::iterator                 Value_iterator;
typedef vector< Vector >                 Vector_cont;
typedef Cartesian_matrix<
  plus< int >,
  Value_iterator,
  Value_iterator >                       Matrix;

int main() {
  // set of vectors the matrices are build from:
  Vector_cont vectors;

  // generate a random vector and sort it:
  Vector a;
  int i;
  cout << "a = ( ";
  for ( i = 0; i < 5; ++i) {
    a.push_back( default_random( 100));
    cout << a.back() << " ";
  }
  cout << ")" << endl;
  sort( a.begin(), a.end(), less< Value >());

  // build a cartesian from a:
  Matrix M( a.begin(), a.end(), a.begin(), a.end());

  // search an upper bound for max(a):
  Value bound( a[4]);
  Value upper_bound(
    sorted_matrix_search(
      &M,
      &M + 1,
      sorted_matrix_search_traits_adaptor(
        bind2nd( greater_equal< Value >(), bound),
        M)));
  cout << "upper bound for " << bound << " is "
       << upper_bound << endl;

  return 0;
} 
