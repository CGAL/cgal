#define CGAL_NEF3_VISUALIZOR

#include <CGAL/basic.h>

#include <CGAL/Gmpz.h>
#include <CGAL/Simple_homogeneous.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_iteration.h>


#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Gmpz                         NT;
typedef CGAL::Simple_homogeneous<NT>       K;
typedef CGAL::SNC_items<K, bool>           SNC_items;
typedef CGAL::SNC_structure<SNC_items>     SNC_structure;
typedef CGAL::Nef_polyhedron_3<SNC_items>  Nef_polyhedron;
typedef CGAL::Polyhedron_3<K>              Polyhedron;

int main() {
    
    Polyhedron P;
    std::cin >> P;
    
    Nef_polyhedron NP(P);
    
    NP.dump();
    NP.visualize();
    
    return 0;
}
