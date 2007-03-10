#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <fstream>

typedef CGAL::Gmpz  NT;
typedef CGAL::Homogeneous<NT>  SK;
typedef CGAL::Extended_homogeneous<NT>  EK;
typedef CGAL::Nef_polyhedron_3<SK>  Nef_polyhedron_S;
typedef CGAL::Nef_polyhedron_3<EK>  Nef_polyhedron_E;

int main() {
  Nef_polyhedron_E E;
  Nef_polyhedron_S S;

  std::cin >> E;

  if(E.is_bounded()) {
    std::ofstream out("temp.nef3");
    out << E;
    std::ifstream in("temp.nef3");
    in >> S;
  }
}
