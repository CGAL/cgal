#include <CGAL/Simple_cartesian.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpq.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Cartesian<double> K2;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> K3;
typedef K3::RT FT3;

int main() {
  K::Point_2 p(-1./3, 2.);
  K::Vector_2 v = p - CGAL::ORIGIN;
  K::Circle_2 c(p, 10);

  K2::Point_2 p2(-1./3, 2.);
  K2::Vector_2 v2 = p2 - CGAL::ORIGIN;
  K2::Circle_2 c2(p2, 10);

  // no correct pretty-printer for CGAL::Gmpq
  K3::Point_3 p3(-3, 10, 2);
  K3::Vector_3 v3 = p3 - CGAL::ORIGIN; 
  return 0;
}
