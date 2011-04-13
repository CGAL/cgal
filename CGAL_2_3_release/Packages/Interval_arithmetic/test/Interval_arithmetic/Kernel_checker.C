
// Test program for the kernel checker.

// Some short names
#define Cartesian CA
#define Cartesian_converter CC
#define Kernel_checker KC
#define Filtered_kernel Fk

#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Kernel_checker.h>

#include <CGAL/Random.h>

#include <CGAL/Delaunay_triangulation_3.h>

typedef CGAL::Filtered_kernel<CGAL::Cartesian<double> > K1;
typedef CGAL::Cartesian<CGAL::MP_Float> K2;
typedef CGAL::Kernel_checker<K1, K2, CGAL::Cartesian_converter<K1, K2> > K;

struct Rep : public  K {};
typedef Rep::RT  NT;
typedef CGAL::Delaunay_triangulation_3<Rep> Delaunay3d;

int my_rand()
{
  return int(CGAL::default_random.get_double()*(1<<31));
}

int main()
{
  Delaunay3d D;
  for (int i=0; i<100; i++)
    D.insert(Rep::Point_3(NT(my_rand()), NT(my_rand()), NT(my_rand())));
  return 0;
}
