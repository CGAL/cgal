

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>


typedef CORE::Expr                                              					NT;					
typedef CGAL::Cartesian<NT>                                     					Kernel;
typedef Kernel::FT                                                                  FT;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> 		Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>				Triangulation;
typedef Triangulation::Face_iterator                                                Face_iterator;

int main(void) {

    Triangulation tr;    
    tr.insert_dummy_points();
    
    cout << "Triangulation successfully initialized with dummy points!" << endl << "---------------------------------------------" << endl;
    cout << "Number of vertices:                  " << tr.number_of_vertices() << endl;
    cout << "Number of faces:                     " << tr.number_of_faces() << endl;
    cout << "Number of edges:                     " << tr.number_of_edges() << endl;
    cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << endl;

    assert(tr.is_valid(true));

    for (Face_iterator fi = tr.faces_begin(); fi != tr.faces_end(); fi++) {
        cout << "Face " << fi->get_number() << " is canonical: " << (fi->is_canonical() ? "TRUE" : "FALSE") << endl;
    }

    return 0;
}