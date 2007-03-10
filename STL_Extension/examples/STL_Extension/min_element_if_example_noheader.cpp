#include <CGAL/algorithm.h>
#include <CGAL/function_objects.h>
#include <vector>
#include <iostream>
#include <functional>

using std::vector;
using std::cout;
using std::endl;
using std::modulus;
using std::greater;
using std::bind2nd;
using CGAL::compose1_1;
using CGAL::min_element_if;

int main()
{
  vector< int > v;
  v.push_back(3);
  v.push_back(5);
  v.push_back(2);
  cout << "min_odd = "
       << *min_element_if(v.begin(),
			  v.end(),
			  compose1_1(bind2nd(greater< int >(), 0),
				     bind2nd(modulus< int >(), 2)))
       << endl;
  return 0;
}
