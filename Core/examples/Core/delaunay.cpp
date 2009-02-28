
#include <CGAL/CORE_Expr.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CORE::Expr Real;
typedef CGAL::Simple_cartesian<Real> K;
typedef CGAL::Delaunay_triangulation_2<K> DT;

typedef K::Point_2 Point_2;

int main() {

  DT dt;
  double two = 2;
  Point_2 p(0,0), q(std::sqrt(two),1), r(0,1);

  dt.insert(p);
  dt.insert(q);
  dt.insert(r);

  std::cout << dt << std::endl;

  return 0;
}


