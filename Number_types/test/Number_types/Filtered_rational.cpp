#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Filtered_rational.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <fstream>

typedef CGAL::Simple_cartesian<CGAL::Filtered_rational> SC;
typedef CGAL::internal::Static_filters<SC> K;
typedef K::Point_2 Point_2;

typedef CGAL::Delaunay_triangulation_2<K>  Triangulation;

int main()
{
  Point_2 p(7.8, 12.);
  Triangulation t;
  t.insert(p);

  CGAL::is_negative(p.x());
  return 0;
}
