// This test a bug in the implementation of the swap of boost variant.
// Here is the bug report:
// https://svn.boost.org/trac/boost/ticket/2839

#include <vector>
#include <boost/variant.hpp>

int main()
{
  typedef boost::variant< std::vector<int>::iterator > Variant;
  std::vector<int> vect;
  Variant x = vect.begin(), y=vect.end();
  x.swap(y);
  return 0;
}
