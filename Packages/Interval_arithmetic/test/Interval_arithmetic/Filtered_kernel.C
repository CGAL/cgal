// Test program for the new filtering wrapper for kernel traits.
// Currently, it's mostly a compilation test.

// We need some short names :
#define Cartesian CA
#define Homogeneous HO
#define Homogeneous_converter HC
#define Cartesian_converter CC
#define Simple_cartesian SC
#define Filtered_kernel Fk

#include <CGAL/Cartesian.h>
// #include <CGAL/Homogeneous.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Homogeneous_converter.h>
#include <CGAL/Cartesian_converter.h>

#include <CGAL/Random.h>
#include <CGAL/MP_Float.h>

#include <CGAL/Delaunay_triangulation_3.h>

#if 1
typedef CGAL::Filtered_kernel<CGAL::Cartesian<double> > Rep;
#else
typedef CGAL::Homogeneous<int> K1;
typedef CGAL::Homogeneous<CGAL::MP_Float> K2;
typedef CGAL::Homogeneous<CGAL::Interval_nt_advanced> K3;

typedef CGAL::Filtered_kernel<K1, K2, K3,
			      CGAL::Homogeneous_converter<K1, K2>,
			      CGAL::Homogeneous_converter<K1, K3> > Rep;
#endif

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
    D.insert(Rep::Point_3(NT(my_rand()), NT(my_rand()), NT(my_rand()),
                          NT(my_rand()) // for homogeneous
			));
  return 0;
}
