#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>


int main()
{
  std::ifstream in("data/double.txt");
  CGAL::Simple_cartesian<double>::Point_3 p;
  while(in >> p){
    std::cout << p << std::endl;
  }
  return 0;
}
