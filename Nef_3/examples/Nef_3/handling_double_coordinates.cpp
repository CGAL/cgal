// file: examples/Nef_3/handling_double_coordinates.C

#include <CGAL/Gmpq.h>
#include <CGAL/Lazy_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

typedef CGAL::Gmpq  NT;
typedef CGAL::Lazy_kernel<CGAL::Simple_cartesian<NT> > Kernel;
typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
typedef CGAL::SNC_indexed_items Items;
typedef CGAL::Nef_polyhedron_3<Kernel, Items>  Nef_polyhedron;

int main() {

  Polyhedron P;
  std::cin >> P;
  Nef_polyhedron N(P);

  std::cout << "Lazy_kernel + Simple_cartesian + SNC_indexed_items"
	    << std::endl
	    << "  allows efficient handling of input "
		 "using floating point coordinates"
	    << std::endl;

  if(N1.is_simple()) {
    N1.convert_to_Polyhedron(P);
    std::cout << P;
  } 
  else {
    std::cout << N1;
  }
}
