#include <iostream>
#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Visibility_complex_2.h>
#include <CGAL/Visibility_complex_2/Circle_traits.h>
#include <CGAL/Gmpz.h>

typedef CGAL::Simple_cartesian<CGAL::Gmpz>          K;
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;
typedef CGAL::Visibility_complex_2<Gt>              VC;
typedef VC::Constraint_input Constraint_input;

int main(int argc, char** argv) {
  std::ifstream di("data/circle.d");
  std::istream_iterator<Gt::Disk> disk_it(di),disk_end;
  std::ifstream ci("data/circle.c");
  std::istream_iterator<Constraint_input> constraint_it(ci),constraint_end;
  VC V(disk_it,disk_end,constraint_it,constraint_end);
  return 0;
}
