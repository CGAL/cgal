// small example for compilation
// check of Delaunay_triangulation_2 using the kernel concept
// archetype

#include <CGAL/basic.h>
#include <CGAL/Triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Kernel_archetype.h>

typedef CGAL::Kernel_archetype  K;
typedef K::Point_2              Point;

typedef CGAL::Delaunay_triangulation_2<K>  Delaunay_triang_2;

Delaunay_triang_2 dt;

int main()
{
  std::list<Point> input;
  Point act;

  dt.insert(act);
  return 0;
}
