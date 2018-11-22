#include <iostream>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Point_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <vector>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>

using namespace CGAL;
using namespace std;

int main(void) {

    typedef CORE::Expr                                                          NT;
    typedef CGAL::Cartesian<NT>                                                 Kernel;
    typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                          CGAL::Hyperbolic_octagon_translation> Traits;
    typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>        Triangulation;

    Triangulation tr;
    CGAL_assertion(tr.is_valid());

    cout << "triangulation works!" << std::endl;
    cout << "nb of vertices: " << tr.number_of_vertices() << endl;
    cout << "nb of faces: " << tr.number_of_faces() << endl;

    return 0;
}