// Computes the minsphere of some random spheres.
// This example illustrates how to use CGAL::Point_3 and CGAL::
// Weighted_point with the Min_sphere_of_spheres_d package.

#include <CGAL/Cartesian.h>
#include <CGAL/Random.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <vector>

const int N = 1000;                       // number of spheres
const int LOW = 0, HIGH = 10000;          // range of coordinates and radii

typedef CGAL::Gmpq                        FT;
//typedef double                          FT;
typedef CGAL::Cartesian<FT>               K;
typedef CGAL::Min_sphere_of_spheres_d_traits_3<K,FT> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
typedef K::Point_3                        Point;
typedef Traits::Sphere                    Sphere;

int main () {
  std::vector<Sphere> S;                  // n spheres
  CGAL::Random r;                         // random number generator

  for (int i=0; i<N; ++i) {
    const FT x = r.get_int(LOW,HIGH),
             y = r.get_int(LOW,HIGH),
             z = r.get_int(LOW,HIGH);
    Point p(x,y,z);                       // random center...
    S.push_back(Sphere(p,r.get_int(LOW,HIGH))); // ...and random radius
  }

  Min_sphere ms(S.begin(),S.end());       // check in the spheres
  CGAL_assertion(ms.is_valid());
}
