#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron;
typedef Nef_polyhedron::Sphere_circle Sphere_circle;

int main()
{
  Nef_polyhedron N1(Sphere_circle(1,0,0));
  Nef_polyhedron N2(Sphere_circle(0,1,0), Nef_polyhedron::EXCLUDED);
  Nef_polyhedron N3 = N1 * N2;
  return 0;
}
