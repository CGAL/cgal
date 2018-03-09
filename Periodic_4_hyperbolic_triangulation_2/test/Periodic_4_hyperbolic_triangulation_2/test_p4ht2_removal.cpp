

#include <boost/tuple/tuple.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_smallint.hpp>
#include <boost/random/variate_generator.hpp>
#include <CGAL/point_generators_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h>
#include <CGAL/Hyperbolic_octagon_translation.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>
#include <CGAL/determinant.h>

typedef CORE::Expr                                                                  NT;
typedef CGAL::Cartesian<NT>                                                         Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel,
                                      CGAL::Hyperbolic_octagon_translation>         Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>                Triangulation;
typedef Triangulation::Face_iterator                                                Face_iterator;
typedef Triangulation::Vertex_handle 												Vertex_handle;
typedef Triangulation::Point 														Point;
typedef Traits::Side_of_original_octagon                                            Side_of_original_octagon;
typedef CGAL::Creator_uniform_2<NT, Point >                                         Creator;

using std::cout;
using std::endl;

int main(int argc, char** argv) {

    if (argc < 2) {
        cout << "usage: " << argv[0] << " [number of vertices to insert]" << endl;
        return -1;
    }

    int N = atoi(argv[1]);

    Triangulation tr;    
    
    cout << "Triangulation successfully initialized with dummy points!" << endl << "---------------------------------------------" << endl;
    cout << "Number of vertices:                  " << tr.number_of_vertices()  << endl;
    cout << "Number of faces:                     " << tr.number_of_faces()     << endl;
    cout << "Number of edges:                     " << tr.number_of_edges()     << endl;
    cout << "Expected edges (by Euler relation):  " << tr.number_of_vertices() + tr.number_of_faces() + 2 << endl << endl;

    assert(tr.is_valid(true));


    CGAL::Random_points_in_disc_2<Point, Creator> g( 1.0 );
    Side_of_original_octagon pred;

    cout << "Inserting " << N << " random new vertices... " << endl;
    std::vector<Vertex_handle> new_v;
    int cnt = 0;
    do {
        Point pt = *g;
        ++g;
        if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
            new_v.push_back(tr.insert(pt));
            cout << "Triangulation has " << tr.number_of_vertices() << " vertices" << endl;
            cnt++;
        }
    } while (cnt < N);
    cout << "Done!" << endl << endl;

    assert(tr.is_valid());
    cout << "Now the triangulation has " << tr.number_of_vertices() << " vertices." << endl << endl;
    cout << "Removing all newly inserted  vertices!" << endl; 

    for (int i = 0; i < new_v.size(); i++) {
        cout << "Removing new vertex " << i << "..." << endl;
        tr.remove(new_v[i]);
    }

    assert(tr.is_valid());
    cout << endl;
    cout << "Now the triangulation has " << tr.number_of_vertices() << " vertices" << endl;

    return 0;
}