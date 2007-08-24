#include <iostream>
#include <fstream>
#include <list>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_2.h>
#include <CGAL/Visibility_complex_2/Circle_traits.h>
#include <CGAL/Gmpz.h>

typedef CGAL::Simple_cartesian<CGAL::Gmpz>          K;
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;
typedef CGAL::Visibility_complex_2<Gt>              VC;

int main()
{
  std::ifstream di("data/circle.d");
  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  // Computing the Visibibility Complex
  VC V(disk_it,disk_end);
  // Printing the bitangents incident to the first disk
  VC::Edge_const_handle e=V.positive_edge(*V.disks_begin());
  VC::Edge_const_handle e1=e;
  do {
    VC::Vertex_const_handle v=e->sup();
    std::cout<<*v<<"\n";
    e=v->ccw_edge(e->object());
  } while (e!=e1);
  return 0;
}
