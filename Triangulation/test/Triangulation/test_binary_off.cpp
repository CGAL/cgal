#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/IO/Triangulation_off_ostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Triangulation;
typedef K::Point_3 Point;

int main()
{
  Triangulation T;
  T.insert(Point(0,0,0));
  T.insert(Point(1,0,0));
  T.insert(Point(0,1,0));
  T.insert(Point(0,0,1));

  const char* filename = "test_binary_output.off";
  std::ofstream out(filename, std::ios::binary);

  if(!out) {
      std::cerr << "Error opening file!" << std::endl;
      return 1;
  }

  out << T;

  out.close();

  return 0;
}
