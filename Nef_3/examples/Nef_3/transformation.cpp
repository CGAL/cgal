#include <CGAL/Exact_integer.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <cassert>


typedef CGAL::Extended_homogeneous<CGAL::Exact_integer>  Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel>  Nef_polyhedron;
typedef Nef_polyhedron::Plane_3  Plane_3;
typedef Nef_polyhedron::Vector_3  Vector_3;
typedef Nef_polyhedron::Aff_transformation_3  Aff_transformation_3;

int main() {

  Nef_polyhedron N(Plane_3(0,1,0,0));
  Aff_transformation_3 transl(CGAL::TRANSLATION, Vector_3(5, 7, 9));
  Aff_transformation_3 rotx90(1,0,0,
                              0,0,-1,
                              0,1,0,
                              1);
  Aff_transformation_3 scale(CGAL::SCALING, 3, 2);

  N.transform(transl);
  assert(N == Nef_polyhedron(Plane_3(0,1,0,-7)));
  N.transform(rotx90);
  assert(N == Nef_polyhedron(Plane_3(0,0,1,-7)));
  N.transform(scale);
  assert(N == Nef_polyhedron(Plane_3(0,0,2,-21)));

  return 0;
}
