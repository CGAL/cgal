#include <CGAL/Exact_integer.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>

typedef CGAL::Homogeneous<CGAL::Exact_integer>  Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3  Vector_3;
typedef Kernel::Aff_transformation_3  Aff_transformation_3;

int main() {
  Polyhedron P;
  std::cin >> P;
  if(P.is_closed()) {
    Nef_polyhedron N1(P);
    Nef_polyhedron N2(N1);
    Aff_transformation_3 aff(CGAL::TRANSLATION, Vector_3(2,2,0,1));
    N2.transform(aff);
    N1 += N2;

    if(N1.is_simple()) {
      N1.convert_to_polyhedron(P);
      std::cout << P;
    }
    else
      std::cerr << "N1 is not a 2-manifold." << std::endl;
  }
}
