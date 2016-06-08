

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
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> 		Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>				Triangulation;

typedef Triangulation::Point_2	                                   					Point_2;
typedef Triangulation::Face_handle													Face_handle;
typedef Triangulation::Vertex_handle												Vertex_handle;
typedef Triangulation::Locate_type 													Locate_type;

typedef std::pair<Vertex_handle, Vertex_handle>                                     Edge;
typedef unsigned short int                                                          Int;
typedef CGAL::Hyperbolic_word_4<Int, Traits>                                        Offset;
typedef std::pair<Edge, bool>                                                       OEdge;
typedef std::set<OEdge>                                                             OEdgeSet;
typedef OEdgeSet::iterator                                                          OEI;

int ccw(int i) {
    return (i+1)%3;
}

int main(void) {

    Triangulation tr;    
    tr.insert_dummy_points();
    
    cout << "Triangulation successfully initialized with dummy points!" << endl;
    cout << "Number of vertices:                  " << tr.number_of_vertices() << endl;
    cout << "Number of faces:                     " << tr.number_of_faces() << endl;
    cout << "Number of edges:                     " << tr.number_of_edges() << endl;
    cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << endl;

    cout << "Triangulation is valid: " << (tr.is_valid(true) ? "YES" : "NO") << endl;

    assert(tr.is_valid());

    return 0;
}