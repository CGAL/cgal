// Computes the minsphere of some random spheres.
#include <CGAL/Cartesian_d.h>
#include <CGAL/Random.h>
#include <CGAL/Min_sphere_of_spheres_d.h>
#include <CGAL/Gmpq.h>
#include <iostream>

typedef CGAL::Gmpq                        FT;
//typedef double                          FT;
typedef CGAL::Cartesian_d<FT>             K;
typedef CGAL::Min_sphere_of_spheres_d_traits_d<K,FT> Traits;
typedef CGAL::Min_sphere_of_spheres_d<Traits> Min_sphere;
typedef K::Point_d                        Point;
typedef Traits::Sphere                    Sphere;

const int n = 1000;                       // number of spheres
const int d = 3;                          // dimension of points
const int Low = 0, High = 10000;          // range of coordinates and radii

int main () {
  Sphere        S[n];                     // n spheres
  FT        coord[d];                     // d coordinates
  CGAL::Random  r;                        // random number generator
  
  for (int i=0; i<n; ++i) {
    for (int j=0; j<d; ++j)
      coord[j] = r.get_int(Low,High);
    Point p(d,coord,coord+d);             // random center...
    S[i] = Sphere(p,r.get_int(Low,High)); // ...and random radius
  }
  
  Min_sphere  ms(S,S+n);                  // smallest enclosing sphere
  CGAL_assertion(ms.is_valid(true,0));
}
