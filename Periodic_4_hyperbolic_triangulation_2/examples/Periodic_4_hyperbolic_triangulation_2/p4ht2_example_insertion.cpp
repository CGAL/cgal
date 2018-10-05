
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
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Triangulation::Point                                                    Point;
typedef Triangulation::size_type                                                size_type;
typedef CGAL::Creator_uniform_2<NT,Point>                                       Creator;


using std::cout;
using std::endl;

int main(int argc, char** argv) {    

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number_of_points_to_insert]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);   

    int N1 = N/2;
    int N2 = N - N1;
    
    // Generate N random points uniformly distributed with respect to the Euclidean 
    // metric in the disk with radius 0.85, which contains the fundamental domain. 
    // Some of the points will be outside the octagon, so they will not be inserted.
    std::vector<Point> pts;
    CGAL::Random_points_in_disc_2<Point,Creator> g( 0.85 );
    CGAL::cpp11::copy_n( g, N1, std::back_inserter(pts) );

    // The triangulation is automatically initialized with the dummy points.
    Triangulation tr;

    // Batch-insert new points in the triangulation. Note that by default, after
    // the insertion of a set of points, an attempt is made to remove all dummy points
    // from the triangulation. Note the third boolean parameter in the call of the 
    // `insert()` function which suppresses the automatic removal of dummy points.
    cout << "Batch-inserting " << N1 << " random points in the triangulation... "; cout.flush();
    size_type N_batch_inserted = tr.insert(pts.begin(), pts.end(), false);
    cout << "DONE! " << endl;

    // Insert new points in the triangulation one by one. When points are inserted 
    // one by one, dummy points are not automatically removed.
    cout << "Single-inserting " << N2 << " random points in the triangulation... "; cout.flush();
    int N_single_inserted = 0;
    for (int i = 0; i < N2; i++) {
        Vertex_handle vh = tr.insert(*(++g));
        if (vh != Vertex_handle()) {
            N_single_inserted++;
        }
    }
    cout << "DONE! " << endl;

    // Total number of inserted points.
    int N_inserted = N_batch_inserted + N_single_inserted;

    // Finally, we try to manually remove all dummy points from the triangulation.
    cout << "Cleaning dummy points from the triangulation... "; cout.flush();
    int DP_remaining = tr.clean_dummy_points();
    cout << "DONE! " << endl;

    // Make sure that the triangulation is valid.
    CGAL_assertion(tr.is_valid());

    int NV = tr.number_of_vertices();
    int NF = tr.number_of_faces();
    int NE = tr.number_of_edges();

    // This function `tr.number_of_dummy_points()` returns the number of dummy points that 
    // are currently in the triangulation. 
    int NDP = tr.number_of_dummy_points();

    // The result must be identical.
    CGAL_assertion(NDP == DP_remaining);

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
