#include <CGAL/Exact_integer.h>
#include <CGAL/Filtered_extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_2.h>

typedef CGAL::Exact_integer RT;
typedef CGAL::Filtered_extended_homogeneous<RT> Extended_kernel;
typedef CGAL::Nef_polyhedron_2<Extended_kernel> Nef_polyhedron;
typedef Nef_polyhedron::Line  Line;

int main()
{
  Nef_polyhedron N1(Line(1,0,0));
  Nef_polyhedron N2(Line(0,1,0), Nef_polyhedron::EXCLUDED);
  Nef_polyhedron N3 = N1 * N2; // line (*)
  return 0;
}
