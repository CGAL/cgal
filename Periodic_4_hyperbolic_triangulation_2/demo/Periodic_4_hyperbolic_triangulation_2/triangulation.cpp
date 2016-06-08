

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>

#include <CGAL/Periodic_4_hyperbolic_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Periodic_4_hyperbolic_triangulation_dummy_14.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>


typedef CORE::Expr                                              					NT;					
typedef CGAL::Cartesian<NT>                                     					Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel> 		Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>				Triangulation;

typedef Triangulation::Point_2	                                   					Point;
typedef Triangulation::Face_handle													Face_handle;
typedef Triangulation::Vertex_handle												Vertex_handle;
typedef Triangulation::Locate_type 													Locate_type;
typedef Traits::FT                                                                  FT;
typedef std::set<Face_handle>::iterator                                             face_iterator;


using namespace CGAL;


int main(void) {

    Triangulation tr;    
    tr.insert_dummy_points();
    std::cout << "Dummy points are locked and loaded!" << std::endl;
    Locate_type lt;
    int li;


    // This is the barycenter of face 10
    Point query = Point(FT(-0.59841), FT(-0.15901));
    Face_handle fh = tr.locate( query, lt, li );
    cout << "Point " << query << " is located in face " << fh->get_number();
    cout << (lt == Triangulation::EDGE ? " [EDGE]" : (lt == Triangulation::VERTEX ? " [VERTEX]" : " [FACE]")) << endl;

    std::set<Face_handle> faces_in_conflict;
    tr.find_in_conflict(query, fh, faces_in_conflict);
    cout << "Faces in conflict: ";
    for (face_iterator it = faces_in_conflict.begin(); it != faces_in_conflict.end(); it++) {
        cout << (*it)->get_number() << ", ";
    }
    cout << endl;

    cout << endl;
    faces_in_conflict.clear();

    // This is the midpoint of edge (0, 1) in face 11
    Point query2(FT(-0.61599), FT(-0.38844));
    fh = tr.locate( query2, lt, li );
    cout << "Point " << query2 << " is located in face " << fh->get_number();
    cout << (lt == Triangulation::EDGE ? " [EDGE]" : (lt == Triangulation::VERTEX ? " [VERTEX]" : " [FACE]")) << endl;
    tr.find_in_conflict(query2, fh, faces_in_conflict);
    cout << "Faces in conflict: ";
    for (face_iterator it = faces_in_conflict.begin(); it != faces_in_conflict.end(); it++) {
        cout << (*it)->get_number() << ", ";
    }
    cout << endl;

    cout << endl;
    faces_in_conflict.clear();

    // This is vertex(0) of face 11
    Point query3(FT(-0.776887), FT(-0.321797));
    fh = tr.locate( query3, lt, li );
    cout << "Point " << query3 << " is located in face " << fh->get_number();
    cout << (lt == Triangulation::EDGE ? " [EDGE]" : (lt == Triangulation::VERTEX ? " [VERTEX]" : " [FACE]")) << endl;
    tr.find_in_conflict(query3, fh, faces_in_conflict);
    cout << "Faces in conflict: ";
    for (face_iterator it = faces_in_conflict.begin(); it != faces_in_conflict.end(); it++) {
        cout << (*it)->get_number() << ", ";
    }
    cout << endl;

    cout << endl;
    faces_in_conflict.clear();

    // This is the origin (duh)
    Point query4(FT(0), FT(0));
    fh = tr.locate( query4, lt, li );
    cout << "Point " << query4 << " is located in face " << fh->get_number();
    cout << (lt == Triangulation::EDGE ? " [EDGE]" : (lt == Triangulation::VERTEX ? " [VERTEX]" : " [FACE]")) << endl;
    tr.find_in_conflict(query4, fh, faces_in_conflict);
    cout << "Faces in conflict: ";
    for (face_iterator it = faces_in_conflict.begin(); it != faces_in_conflict.end(); it++) {
        cout << (*it)->get_number() << ", ";
    }
    cout << endl;    

    cout << endl << "Triangulation is valid: " << (tr.is_valid() ? "TRUE" : "FALSE") << endl;

    return 0;
}