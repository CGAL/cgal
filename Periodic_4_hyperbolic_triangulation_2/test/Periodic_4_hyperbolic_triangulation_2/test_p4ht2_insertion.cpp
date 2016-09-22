

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>

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

typedef Kernel::Point_2                                                             Point;
typedef CGAL::Creator_uniform_2<double, Point>                                      Creator;
typedef std::vector<Point>                                                          Vector;
typedef Triangulation::Face_handle                                                  Face_handle;
typedef Triangulation::Vertex_handle                                                Vertex_handle;
typedef Triangulation::Locate_type                                                  Locate_type;
typedef Triangulation::Offset                                                       Offset;

int ccw(int i) {
    return (i+1)%3;
}

int main(void) {

    Triangulation tr;    
    tr.insert_dummy_points();

    int N = 100;
    Vector pts;
    pts.reserve(N);
    CGAL::Random_points_in_disc_2<Point, Creator> g( 1.0 );
    CGAL::cpp11::copy_n( g, N, std::back_inserter(pts));

    Locate_type lt;
    int li;

    int bad = 0;
    for (int i = 0; i < N; i++) {
        Offset loff;
        Face_handle fh = tr.euclidean_visibility_locate( pts[i], lt, li, loff );
        Vertex_handle vh = tr.insert(pts[i], fh);
        if (vh == Vertex_handle())
            bad++;
    }

    cout << "Tried to insert " << N << " random points, " << bad << " were rejected." << endl;
    cout << "Number of vertices:                  " << tr.number_of_vertices() << endl;
    cout << "Number of faces:                     " << tr.number_of_faces() << endl;
    cout << "Number of edges:                     " << tr.number_of_edges() << endl;
    cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << endl;

    cout << "Triangulation is valid: " << (tr.is_valid(true) ? "YES" : "NO") << endl;

    assert(tr.is_valid());

    

    return 0;
}