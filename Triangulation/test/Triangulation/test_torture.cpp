#include <CGAL/config.h>
#if defined(BOOST_GCC) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)

#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else

#include <CGAL/Epick_d.h>
#include <CGAL/internal/Combination_enumerator.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/algorithm.h>
#include <tilted_grid.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

// Inserts N^D points then removes min(N^D, 100) points
template<typename DC >
void test(const int D, const int d, const int N, bool no_transform)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that DC contains sub-types, before
    // instanciating the function.
    typedef typename DC::Point Point;
    typedef typename DC::Geom_traits::RT RT;

    DC dc(D);

    vector<RT> coords(D);
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
        int nb = rand() % 1000;
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
                    c = (rand() % 11) - 5;
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
                coords[i] = coords[i] + (*pit)[j] * aff[j][i];
        }
        points.push_back(Point(D, coords.begin(), coords.end()));
    }
    assert( dc.is_valid() );
    cout << " Inserting " << points.size() << " points.";
    dc.insert(points.begin(), points.end());
    assert( d >= dc.current_dimension() );
    assert( points.size() >= dc.number_of_vertices() );
    if( points.size() > dc.number_of_vertices() )
        assert( d > dc.current_dimension() );
    assert( dc.is_valid() );
    if( 2 == dc.current_dimension() )
        assert( 2 * dc.number_of_vertices() == dc.number_of_full_cells() + 2 );
    if( dc.current_dimension() > 3 )
    {
        CGAL::cpp98::random_shuffle(points.begin(), points.end());
        if (points.size() > 100)
          points.resize(100);
    }
    cout << " Removing " << points.size() << " points.";
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
    typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > K;
    //typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
    typedef CGAL::Delaunay_triangulation<K> Triangulation;
    for( int d = 0; d <= D; ++d )
    {
        cout << "\n\n** Delaunay of " << d
             << "-dimensional regular grid in R^" << D << endl;
        for( int t = 0; t < nb_trials; ++t )
        {
            cout << "\n Run " << (t + 1) << ": ";
            test<Triangulation>(D, d, N, t < nb_trials / 2);
        }
    }
}

int main(int argc, char **argv)
{
    int N = 3;
    int nb_trials = 2;
    unsigned int rand_init = static_cast<unsigned int>(time(nullptr));
    if( argc > 1 )
        N = atoi(argv[1]);
    if( argc > 2 )
        nb_trials = atoi(argv[2]);
   if( argc > 3 )
        rand_init = atoi(argv[3]);
   std::cout<<argv[0]<<" "<<N<<" "<<nb_trials<<" "<<rand_init<<std::endl;
    srand(rand_init);
    go<1>(N, nb_trials);
    go<2>(N, nb_trials);
    go<3>(N, nb_trials);
    go<4>(N, nb_trials);
#ifdef NDEBUG
    go<5>(N, nb_trials);
#endif
    cout << std::endl;
    return 0;
}

#endif
