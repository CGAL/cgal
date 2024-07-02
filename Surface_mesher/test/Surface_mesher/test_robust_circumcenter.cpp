#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Robust_circumcenter_traits_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K1;
typedef CGAL::Robust_circumcenter_traits_3<K1>  K;
typedef CGAL::Delaunay_triangulation_3<K>   Triangulation;
typedef K::Point_3                          Point;

int main()
{
  Triangulation tr;
  int a, b, d;
  for (a=0;a!=4;a++)
    for (b=0;b!=4;b++)
      for (d=0;d!=4;d++)
        tr.insert(Point((a*b-d*a)*10 +a ,(a-b+d +5*b)*100,
                            a*a-d*d-b));
  Triangulation::Finite_cells_iterator cit=tr.finite_cells_begin();
  for(; cit != tr.finite_cells_end(); ++cit) {
    Point circum = tr.dual(cit);
    CGAL_USE(circum);
  }
  return 0;
}
