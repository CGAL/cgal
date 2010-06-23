#ifndef NEF_POLYHEDRON_S2_CREATE_RANDOM_H
#define NEF_POLYHEDRON_S2_CREATE_RANDOM_H

#include <CGAL/random_selection.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Nef_polyhedron_S2.h>

namespace CGAL {

template <typename K,typename I,typename M>
void
create_random_Nef_S2(Nef_polyhedron_S2<K,I,M>& P, int n=5, int seed=0) {

  typedef Nef_polyhedron_S2<K,I,M> Polyhedron;
  typedef typename Polyhedron::Sphere_circle Sphere_circle;
  typedef typename K::Point_3 Point_3;
  typedef typename K::Plane_3 Plane_3;
  typedef typename K::RT RT;

  typedef CGAL::Creator_uniform_3<RT,Point_3>  Creator;
  typedef CGAL::Random_points_in_cube_3<Point_3,Creator> Point_source;

  if(seed == 0)
    std::srand(static_cast<unsigned int>(time(0)));
  else
    std::srand(seed);

  std::list<Sphere_circle> L;
  Point_source S(5);
  Point_3 ph;
  Point_3 o(0,0,0);
  while ( n-- > 0 ) {
    do { ph = *S++; } 
    while ( ph == o );
    Plane_3 h(o,(ph-CGAL::ORIGIN).direction());
    L.push_back( Sphere_circle(h) );
  }

  P = Polyhedron(L.begin(), L.end(), 0.5);

}

} //namespace CGAL
#endif // NEF_POLYHEDRON_S2_CREATE_RANDOM_H
