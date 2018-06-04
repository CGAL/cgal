

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/internal/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Algebraic_kernel_for_circles_2_2.h>
#include <CGAL/Circular_kernel_2.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>


typedef CORE::Expr                                              				NT;				
typedef CGAL::Algebraic_kernel_for_circles_2_2<NT>                              AK;   // Algebraic kernel  
typedef CGAL::Cartesian<NT>                                                     BK;   // Basic kernel
typedef CGAL::Circular_kernel_2<BK, AK>                                         Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
									  CGAL::Hyperbolic_octagon_translation>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Triangulation::Face_iterator                                            Face_iterator;

using std::cout;
using std::endl;

int main(void) {

    Triangulation tr;    
    
    cout << "Triangulation successfully initialized with dummy points!" << endl << "---------------------------------------------" << endl;
    cout << "Number of vertices:                  " << tr.number_of_vertices() << endl;
    cout << "Number of faces:                     " << tr.number_of_faces() << endl;
    cout << "Number of edges:                     " << tr.number_of_edges() << endl;
    cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << endl;

    assert(tr.is_valid(true));

    return 0;
}