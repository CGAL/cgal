#include <CGAL/point_generators_d.h>
//#include <CGAL/Simple_cartesian_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Delaunay_complex.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/algorithm.h>
#include <CGAL/Gmpq.h>
#include <tilted_grid.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

template<typename DC>
void test(const int d, const string & type, const int N)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that DC contains sub-types, before
    // instanciating the function.
    typedef typename DC::Vertex Vertex;
    typedef typename DC::Vertex_handle Vertex_handle;
    typedef typename DC::Simplex Simplex;
    typedef typename DC::Simplex_handle Simplex_handle;
    typedef typename DC::Facet Facet;
    typedef typename DC::Face Face;
    typedef typename DC::Point Point;
    typedef typename DC::Geom_traits::RT RT;
    typedef typename DC::Finite_simplex_const_iterator Finite_simplex_const_iterator;

    typedef CGAL::Random_points_in_iso_box_d<Point> Random_points_iterator;

    DC pc(d);
    cerr << "\nChecking Delaunay_complex of (" << type << d << ") dimension with " << N << " points";
    assert(pc.empty());

    vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 2.0, rng);
    //CGAL::copy_n(rand_it, N, std::back_inserter(points));
    
	vector<int> coords(d);
    for( int i = 0; i < N; ++i )
    {
        for( int j = 0; j < d; ++j )
            coords[j] = lrand48() % 100000;
        points.push_back(Point(d, coords.begin(), coords.end()));
    }
    pc.insert(points.begin(),  points.end());
    cerr << "\nChecking topology and geometry...";
    assert( pc.is_valid() );

    cerr << "\nTraversing finite simplices... ";
    size_t nbfs(0), nbis(0);
    Finite_simplex_const_iterator fsit = pc.finite_simplices_begin();
    while( fsit != pc.finite_simplices_end() )
        ++fsit, ++nbfs;
    cerr << nbfs << " + ";
    vector<Simplex_handle> infinite_simplices;
    pc.gather_incident_simplices(pc.infinite_vertex(), std::back_inserter(infinite_simplices));
    nbis = infinite_simplices.size();
    cerr << nbis << " = " << (nbis+nbfs)
    << " = " << pc.number_of_simplices();
    cerr << "\nComplex has current dimension " << pc.current_dimension();

    // Count convex hull vertices:
    typedef std::vector<Face> Faces;
    Faces edges;
    std::back_insert_iterator<Faces> out(edges);
    pc.gather_incident_upper_faces(pc.infinite_vertex(), 1, out);
    // Count the number of points on the convex hull
    std::cout << "\nThere are " << edges.size() << " vertices on the convex hull.";
    edges.clear();

    // Remove all !
    std::cerr << "\nBefore removal: " << pc.number_of_vertices() << " vertices. After: ";
    std::random_shuffle(points.begin(),  points.end());
    pc.remove(points.begin(),  points.end());
    assert( pc.is_valid() );
    std::cerr << pc.number_of_vertices() << " vertices.";
    // assert( pc.empty() ); NOT YET !
    // CLEAR
    pc.clear();
    assert( -1 == pc.current_dimension() );
    assert( pc.empty() );
    assert( pc.is_valid() );
}

template< int D >
void go(const int N)
{
    typedef double RT;
    //typedef CGAL::Gmpq RT;
    typedef CGAL::Cartesian_d<RT> K;
    //typedef CGAL::Simple_cartesian_d<RT, D> K;
    typedef CGAL::Filtered_kernel_d<K> FK;
    typedef CGAL::Delaunay_complex<FK> Triangulation;
    test<Triangulation>(D, "dynamic", N);
    //test<Triangulation>(D, "static", N);
}

int main(int argc, char **argv)
{
    srand48(time(NULL));
    int N = 100;
    if( argc > 1 )
        N = atoi(argv[1]);
     go<5>(N);
     go<4>(N);
     go<3>(N);
     go<2>(N);
     go<1>(N);

    cerr << std::endl;
    return 0;
}
