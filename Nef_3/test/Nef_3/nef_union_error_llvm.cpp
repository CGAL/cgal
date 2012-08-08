#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <iostream>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::Polyhedron_3<CGAL::Simple_cartesian<double> > Polyhedron_d;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;

int main(int, char**) {
    Polyhedron P1;
    Polyhedron P2;
    Polyhedron P3;
   
    std::ifstream f1("data/icosahedron.off");
    f1 >> P1;
    std::ifstream f2("data/icosahedron.off");
    f2 >> P2;
    Nef_polyhedron N1(P1);
    Nef_polyhedron N2(P2);

    Aff_transformation_3 aff(CGAL::TRANSLATION, Vector_3(1,1,0,1));
    N2.transform(aff);
   
    std::cout << N1.number_of_volumes() << std::endl;
    std::cout << N2.number_of_volumes() << std::endl;
   
    N1+=N2;  // model-dependent crash here
   
    std::cout << N1.number_of_volumes() << std::endl;   
    N1.closure().convert_to_polyhedron(P3);
   
    return 0;
}
