#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_set_3.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;


int main()
{
  Point_set ps;
  ps.insert(Point(0,0,0));
  ps.insert(Point(1,1,1));

  Point_set::Property_map<unsigned char>
    red   = ps.add_property_map<unsigned char>("red"  , 0).first,
    green = ps.add_property_map<unsigned char>("green", 0).first,
    blue  = ps.add_property_map<unsigned char>("blue" , 0).first;


  int i = 1;
  for (Point_set::iterator it = ps.begin(); it != ps.end(); ++ it){
    put(red,   *it, static_cast<unsigned char>(11 * i++));
    put(green, *it, static_cast<unsigned char>(11 * i++));
    put(blue,  *it, static_cast<unsigned char>(11 * i++));
  }

  {
    std::ofstream out("ascii.ply");
    out << ps;
    out.close();

    Point_set ps2;
    std::ifstream in("ascii.ply");
    in >> ps2;

    std::cout << ps2 << std::endl;
  }

  {
    std::ofstream out("binary.ply", std::ios::binary);
    CGAL::IO::set_binary_mode(out);
    out << ps;
    out.close();

    Point_set ps2;
    std::ifstream in("binary.ply", std::ios::binary);
    CGAL::IO::set_binary_mode(in);
    in >> ps2;

    std::cout << ps2 << std::endl;
  }


  return 0;
}
