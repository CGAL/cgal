#include <CGAL/internal/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>
#define USE_NEW_KERNEL
#ifndef USE_NEW_KERNEL
#include <CGAL/Cartesian_d.h> // this is for Old_kernel_d
#else
#include <CGAL/Simple_cartesian_d.h> // this is for New_kernel_d
#endif
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Delaunay_complex.h>
#include <CGAL/algorithm.h>
#include <tilted_grid.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

template<typename DC >
void test(const int D, const int d, const int N, bool no_transform)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that DC contains sub-types, before
    // instanciating the function.
    typedef typename DC::Vertex Vertex;
    typedef typename DC::Vertex_handle Vertex_handle;
    typedef typename DC::Simplex Simplex;
    typedef typename DC::Simplex_handle Simplex_handle;
    typedef typename DC::Facet Facet;
    typedef typename DC::Point Point;
    typedef typename DC::Geom_traits::RT RT;
    typedef typename DC::Finite_simplex_const_iterator Finite_simplex_const_iterator;

    DC dc(D);

    vector<int> coords(D);
    vector<Point> points;
    CGAL::Random rng;

    typedef Tilted_grid<typename DC::Geom_traits> Grid;
    Grid grid;
    typedef typename Grid::PVec PVec;
    PVec pvec;
    grid(d, N, pvec);
    points.clear();

    int k = (d == 0) ? 1 : d;
    CGAL::internal::Combination_enumerator combi(k, 0, D-1);
    if( no_transform )
    {   // choose a random set of axes:
        int nb = lrand48() % 1000;
        for( int i = 0; i < nb; ++i )
        {
            ++combi;
            if( combi.end() )
                combi.init();
        }
    }
    vector<vector<int> > aff;
    for( int j = 0; j < d; ++j )
    {
        aff.push_back(vector<int>());
        for( int i = 0; i < D; ++i )
        {
            if( no_transform )
                aff[j].push_back((combi[j]==i)?1:0);
            else
            {
                int c(0);
                while( 0 == c )
                    c = (lrand48() % 11) - 5;
                aff[j].push_back(c);
            }
        }
    }
    for( typename PVec::iterator pit = pvec.begin(); pit != pvec.end(); ++pit )
    {
        for( int i = 0; i < D; ++i )
        {
            coords[i] = 0;
            for( int j = 0; j < d; ++j )
                coords[i] += (*pit)[j] * aff[j][i];
        }
#ifdef USE_NEW_KERNEL
        points.push_back(Point(coords)); // this is for New_kernel_d
#else
        points.push_back(Point(D, coords.begin(), coords.end())); // this is for Old_kernel_d
#endif
    }
    assert( dc.is_valid() );
    dc.insert(points.begin(), points.end());
    assert( d >= dc.current_dimension() );
    assert( points.size() >= dc.number_of_vertices() );
    if( points.size() > dc.number_of_vertices() )
        assert( d > dc.current_dimension() );
    assert( dc.is_valid() );
    if( 2 == dc.current_dimension() )
        assert( 2 * dc.number_of_vertices() == dc.number_of_simplices() + 2 );
    if( dc.current_dimension() > 4 )
    {
        std::random_shuffle(points.begin(), points.end());
        points.resize(100);
    }
    dc.remove(points.begin(), points.end());
    assert( dc.is_valid() );
    dc.clear();
    assert( -1 == dc.current_dimension() );
    assert( dc.empty() );
    assert( dc.is_valid() );
}

template< int D >
void go(const int N, const int nb_trials)
{
    typedef double RT;
#ifndef USE_NEW_KERNEL
    typedef CGAL::Cartesian_d<RT> K; // this is for Old_kernel_d
#else
    typedef CGAL::Simple_cartesian_d<RT, D> K; // this is for New_kernel_d
#endif
    typedef CGAL::Filtered_kernel_d<K> FK;
    typedef CGAL::Delaunay_complex<FK> Triangulation;
    for( int d = 0; d <= D; ++d )
    {
        cerr << "\nDelaunay of " << d
             << "-dimensional regular grid in R^" << D;
        for( int t = 0; t < nb_trials; ++t )
        {
            cerr << "\n Run " << (t + 1) << ": ";
            test<Triangulation>(D, d, N, t < nb_trials / 2);
        }
    }
}

int main(int argc, char **argv)
{
    srand48(time(NULL));
    int N = 5;
	int nb_trials = 4;
    if( argc > 1 )
        N = atoi(argv[1]);
    if( argc > 2 )
        nb_trials = atoi(argv[2]);
    //go<1>(N, nb_trials);
    //go<2>(N, nb_trials);
    //go<3>(N, nb_trials);
    go<4>(N, nb_trials);
    go<5>(N, nb_trials);
    cerr << std::endl;
    return 0;
}
