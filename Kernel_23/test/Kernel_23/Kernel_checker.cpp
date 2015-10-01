// Test program for the kernel checker.

#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Kernel_checker.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Random.h>

#include <CGAL/Delaunay_triangulation_3.h>

// Does not fully work yet, so I choose a simple case K1 == K2 :)

//typedef CGAL::Filtered_kernel<CGAL::Cartesian<double> > K1;
typedef CGAL::Cartesian<CGAL::MP_Float> K2;
typedef K2 K1;
// typedef CGAL::Cartesian<int> K2;

//typedef CGAL::Kernel_checker<K1, K2, CGAL::Cartesian_converter<K1, K2> > K;
typedef CGAL::Kernel_checker<K1, K2> K;

typedef K::RT  NT;
typedef CGAL::Delaunay_triangulation_3<K> Delaunay3d;

int my_rand()
{
  return int(CGAL::get_default_random().get_double()*(1<<31));
}

int main()
{
  Delaunay3d D;
  for (int i=0; i<100; i++) {
    double x = my_rand(), y = my_rand(), z = my_rand();
    D.insert(K::Point_3(K1::Point_3(x, y, z), K2::Point_3(x, y, z)));
  }
  return 0;
}
