// examples/Nef_S2/visualization.C
// -------------------------------------
#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Nef_polyhedron_S2.h>
#include <CGAL/IO/Nef_polyhedron_S2_OGLUT_stream.h>
#include <CGAL/Nef_S2/create_random_Nef_S2.h>

typedef CGAL::Gmpz RT;
typedef CGAL::Homogeneous<RT> Kernel;
typedef CGAL::Nef_polyhedron_S2<Kernel> Nef_polyhedron_S2;

int main() {

  int n;
  std::cin >> n;

  Nef_polyhedron_S2 S;
  create_random_Nef_S2(S,n);

  CGAL::ogl << S;
  CGAL::ogl << "Visualization example"; 
  CGAL::ogl.display();

  return 0;
}


