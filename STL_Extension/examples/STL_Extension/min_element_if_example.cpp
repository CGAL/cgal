#include <CGAL/algorithm.h>
#include <CGAL/function_objects.h>
#include <vector>
#include <iostream>
#include <boost/functional.hpp>


int main()
{
  std::vector< int > v;
  v.push_back(3);
  v.push_back(5);
  v.push_back(2);
  std::cout << "min_odd = "
            << *CGAL::min_element_if(v.begin(),
                                     v.end(),
                                     [](int i){ return (i%2) > 0; })
            << std::endl;
  return 0;
}
