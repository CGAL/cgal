
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/point_generators_2.h>

typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Kernel::Point_2                                                         Point;
typedef CGAL::Creator_uniform_2<NT,Point>                                       Creator;
typedef Triangulation::size_type                                                size_type;

using std::cout;
using std::endl;

int main(int argc, char** argv) {    

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_points_to_insert]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);   
    
    // Generate N random points uniformly distributed with respect to the Euclidean 
    // metric in the disk with radius 0.85, which contains the fundamental domain. 
    // Some of the points will be outside the octagon, so they will not be inserted.
    std::vector<Point> pts;
    CGAL::Random_points_in_disc_2<Point,Creator> g( 0.85 );
    CGAL::cpp11::copy_n( g, N, std::back_inserter(pts) );

    // The triangulation is automatically initialized with the dummy points.
    Triangulation tr;

    // Insert the new points in the triangulation.
    cout << "Inserting " << N << " random points in the triangulation... "; cout.flush();
    size_type N_inserted = tr.insert(pts.begin(), pts.end());
    cout << "DONE! " << endl;

    // Make sure that the triangulation is valid.
    CGAL_assertion(tr.is_valid());

    int NV = tr.number_of_vertices();
    int NF = tr.number_of_faces();
    int NE = tr.number_of_edges();

    // This function `tr.number_of_dummy_points()` returns the number of dummy points that 
    // are currently in the triangulation. 
    int NDP = tr.number_of_dummy_points();

    cout << endl;
    cout << "-------------- STATS --------------" << endl;
    cout << "Random points generated:           " << N << endl;
    cout << "Vertices in the triangulation:     " << NV << endl;
    cout << "Dummy points in the triangulation: " << NDP << endl;
    cout << "Random points inserted:            " << N_inserted << endl;
    cout << "Random points outside/duplicates:  " << (N - N_inserted) << endl;
    cout << endl;
    cout << "---------- COMBINATORICS ----------" << endl;
    cout << "Number of vertices NV:             " << NV << endl;
    cout << "Number of faces NF:                " << NF << endl;
    cout << "Number of edges NE:                " << NE << endl;

    // The number of vertices in the triangulation must equal the number of points
    // inserted, plus the dummy points in the triangulation.
    CGAL_assertion( N_inserted + NDP == NV );
    cout << "Number of vertices is correct!     " << endl;    

    // Note that the Euler relation is already verified by the function `is_valid()`.
    CGAL_assertion( (2 + NE - NV - NF) / 2 == 2);
    cout << "Euler relation verified!           " << endl << endl;    

    return 0;
}
