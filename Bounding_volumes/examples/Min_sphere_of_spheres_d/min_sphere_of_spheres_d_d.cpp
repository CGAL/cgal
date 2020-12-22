// Computes the minsphere of some random spheres.

#include <CGAL/Cartesian_d.h>
#include <CGAL/Random.h>
#include <CGAL/Exact_rational.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <vector>

const int N = 1000;                       // number of spheres
const int D = 3;                          // dimension of points
const int LOW = 0, HIGH = 10000;          // range of coordinates and radii

typedef CGAL::Exact_rational              FT;
//typedef double                          FT;
typedef CGAL::Cartesian_d<FT>             K;
typedef CGAL::Min_sphere_of_spheres_d_traits_d<K,FT,D> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
typedef K::Point_d                        Point;
typedef Traits::Sphere                    Sphere;

int main () {
  std::vector<Sphere> S;                  // n spheres
  FT coord[D];                            // d coordinates
  CGAL::Random r;                         // random number generator

  for (int i=0; i<N; ++i) {
    for (int j=0; j<D; ++j)
      coord[j] = r.get_int(LOW,HIGH);
    Point p(D,coord,coord+D);             // random center...
    S.push_back(Sphere(p,r.get_int(LOW,HIGH))); // ...and random radius
  }

  Min_sphere ms(S.begin(),S.end());       // check in the spheres
  CGAL_assertion(ms.is_valid());
}
