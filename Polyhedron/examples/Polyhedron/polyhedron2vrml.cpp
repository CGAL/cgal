#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_VRML_1_ostream.h>
#include <iostream>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel>     Polyhedron;

int main() {
    Polyhedron P;
    std::cin >> P;
    CGAL::VRML_1_ostream out( std::cout);
    out << P;
    return ( std::cin && std::cout) ? 0 : 1;
}
