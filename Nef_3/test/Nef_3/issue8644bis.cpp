#include <CGAL/Extended_cartesian.h>
#include <CGAL/Gmpq.h>
#include <CGAL/Nef_polyhedron_3.h>

using Kernel = CGAL::Extended_cartesian<CGAL::Gmpq>;
using Nef = CGAL::Nef_polyhedron_3<Kernel>;

int main()
{
  Nef hspace_1(Nef::Plane_3(1.0, 0.0, 0.0, 0.0), Nef::INCLUDED);
  Nef hspace_2(Nef::Plane_3(1.0, 0.0, 0.0, 1.0), Nef::INCLUDED);
  Nef hspace_3(Nef::Plane_3(0.0, 0.0, 1.0, 1.0), Nef::INCLUDED);

  Nef intersection_1 = hspace_1*hspace_2;
  Nef intersection_2 = hspace_2*hspace_3;

  return 0;
}
