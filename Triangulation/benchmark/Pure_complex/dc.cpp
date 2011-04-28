#include <CGAL/point_generators_d.h>
#include <CGAL/Cartesian_d.h>
//#include <CGAL/Simple_cartesian_d.h>
//#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Delaunay_complex.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

template<typename DC>
void test(const int d, const string & type, const int N)
{
    typedef typename DC::Vertex Vertex;
    typedef typename DC::Vertex_handle Vertex_handle;
    typedef typename DC::Simplex Simplex;
    typedef typename DC::Simplex_handle Simplex_handle;
    typedef typename DC::Facet Facet;
    typedef typename DC::Point Point;
    typedef typename DC::Geom_traits::RT RT;
    typedef typename DC::Finite_simplex_const_iterator Finite_simplex_const_iterator;

    typedef CGAL::Random_points_in_iso_box_d<Point> Random_points_iterator;

    DC pc(d);
    cout << "\nBench'ing Delaunay_complex of (" << type << d << ") dimension with " << N << " points";
    assert(pc.empty());

    vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 2.0, rng);
    CGAL::copy_n(rand_it, N, std::back_inserter(points));
    pc.insert(points.begin(), points.end());
    int nbfs(0);
    cout << '\n' << pc.number_of_vertices() << " vertices, " << (nbfs = pc.number_of_finite_simplices())
        << " finite simplices and " << (pc.number_of_simplices() - nbfs) << " convex hull Facets.";
}

template< int D, typename RT >
void go(const int N)
{
    //typedef CGAL::Simple_cartesian_d<RT, D> K;
    typedef CGAL::Cartesian_d<RT> K;
    //typedef CGAL::Filtered_kernel_d<K> FK;
    typedef CGAL::Delaunay_complex<K/*FK*/> Triangulation;
    test<Triangulation>(D, "static", N);
}

int main(int argc, char **argv)
{
    srand48(time(NULL));
    int N = 10;
    if( argc > 1 )
        N = atoi(argv[1]);
    // go<7>(N);
    // go<6>(N);
    // go<5>(N);
    go<4, double>(N);
    // go<3>(N);
    // go<2>(N);
    // go<2>(N);
    // go<1>(N);

    cout << std::endl;
    return 0;
}
