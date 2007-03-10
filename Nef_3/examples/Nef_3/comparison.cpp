#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>

typedef CGAL::Gmpz  NT;
typedef CGAL::Extended_homogeneous<NT>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
//typedef Nef_polyhedron::Plane_3  Plane_3;
typedef Kernel::Plane_3 Plane_3;
int main() {

  Nef_polyhedron N1(Plane_3(2,5,7,11), Nef_polyhedron::INCLUDED);
  Nef_polyhedron N2(Plane_3(2,5,7,11), Nef_polyhedron::EXCLUDED);

  CGAL_assertion(N1 >= N2);
  CGAL_assertion(N2 <= N1);
  CGAL_assertion(N1 != N2);
  CGAL_assertion(N1 > N2);
  CGAL_assertion(N2 < N1);

  N2 = N2.closure();
  CGAL_assertion(N1==N2);
  CGAL_assertion(N1>=N2);
  CGAL_assertion(N1<=N2);

  return 0;
}
