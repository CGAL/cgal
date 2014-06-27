#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>

typedef CGAL::Exact_integer RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Sphere_point Sphere_point;
typedef Nef_polyhedron::Sphere_segment Sphere_segment;
typedef Nef_polyhedron::Sphere_circle Sphere_circle;

int main() {

  Nef_polyhedron N1(Nef_polyhedron::COMPLETE);

  Sphere_circle c(1,1,1); // c : x + y + z = 0
  Nef_polyhedron N2(c, Nef_polyhedron::INCLUDED);
  Nef_polyhedron N3(N2.complement());
  CGAL_assertion(N1 == N2.join(N3));

  Sphere_point   p1(1,0,0), p2(0,1,0), p3(0,0,1);
  Sphere_segment s1(p1,p2), s2(p2,p3), s3(p3,p1);
  Sphere_segment triangle[3] = { s1, s2, s3 };
  Nef_polyhedron N4(triangle, triangle+3);
  Nef_polyhedron N5;
  N5 += N2;
  N5 = N5.intersection(N4);
  CGAL_assertion(N5 <= N2 && N5 != N4);

  return 0;
}
