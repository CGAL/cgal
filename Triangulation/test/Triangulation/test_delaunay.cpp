#include <CGAL/config.h>
#if defined(BOOST_GCC) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)

#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else

#include <CGAL/Epick_d.h>
#include <CGAL/Epeck_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>
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
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that DC contains sub-types, before
    // instanciating the function.
    typedef typename DC::Full_cell_handle Full_cell_handle;
    typedef typename DC::Face Face;
    typedef typename DC::Point Point;
    typedef typename DC::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;
    typedef typename DC::Finite_vertex_iterator Finite_vertex_iterator;

    DC dt(d);
    cerr << "\nBuilding Delaunay triangulation of (" << type << d << ") dimension with " << N << " points";
    assert(dt.empty());

    vector<Point> points;
    //CGAL::Random rng;
    //Random_points_iterator rand_it(d, 2.0, rng);
    //std::copy_n(rand_it, N, back_inserter(points));

    srand(10);
    for( int i = 0; i < N; ++i )
    {
        vector<double> coords(d);
        for( int j = 0; j < d; ++j )
            coords[j] = static_cast<double>(rand() % 100000)/10000;
        points.push_back(Point(d, coords.begin(), coords.end()));
    }
    dt.insert(points.begin(),  points.end());
    cerr << "\nChecking topology and geometry...";
    assert( dt.is_valid() );

    cerr << "\nTraversing finite full_cells... ";
    size_t nbfs(0), nbis(0);
    Finite_full_cell_const_iterator fsit = dt.finite_full_cells_begin();
    while( fsit != dt.finite_full_cells_end() )
        ++fsit, ++nbfs;
    cerr << nbfs << " + ";
    vector<Full_cell_handle> infinite_full_cells;
    dt.tds().incident_full_cells(dt.infinite_vertex(), back_inserter(infinite_full_cells));
    nbis = infinite_full_cells.size();
    cerr << nbis << " = " << (nbis+nbfs)
    << " = " << dt.number_of_full_cells();
    cerr << "\nThe triangulation has current dimension " << dt.current_dimension();
    CGAL_assertion( dt.number_of_full_cells() == nbis+nbfs);

    cerr << "\nTraversing finite vertices... ";
    size_t nbfv(0);
    Finite_vertex_iterator fvit = dt.finite_vertices_begin();
    while( fvit != dt.finite_vertices_end() )
        ++fvit, ++nbfv;
    cerr << nbfv <<endl;

    // Count convex hull vertices:
    if( dt.maximal_dimension() > 1 )
    {
        typedef vector<Face> Faces;
        Faces edges;
        back_insert_iterator<Faces> out(edges);
        dt.tds().incident_faces(dt.infinite_vertex(), 1, out);
        cout << "\nThere are " << edges.size() << " vertices on the convex hull.";
        edges.clear();
    }
    else // dt.maximal_dimension() == 1
    {
        typedef vector<Full_cell_handle> Cells;
        Cells cells;
        back_insert_iterator<Cells> out(cells);
        dt.tds().incident_full_cells(dt.infinite_vertex(), out);
        cout << "\nThere are " << cells.size() << " vertices on the convex hull.";
        cells.clear();
    }

    // Remove all !
    cerr << "\nBefore removal: " << dt.number_of_vertices() << " vertices. After: ";
    CGAL::cpp98::random_shuffle(points.begin(),  points.end());
    dt.remove(points.begin(),  points.end());
    assert( dt.is_valid() );
    cerr << dt.number_of_vertices() << " vertices.";
    // assert( dt.empty() ); NOT YET !
    // CLEAR
    dt.clear();
    assert( -1 == dt.current_dimension() );
    assert( dt.empty() );
    assert( dt.is_valid() );
}

template< int D >
void go(const int N)
{
    typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > KI;
    typedef CGAL::Delaunay_triangulation<KI> Triangulation;
    test<Triangulation>(D, "inexact static", N);

    typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> KI_dyn;
    typedef CGAL::Delaunay_triangulation<KI_dyn> Triangulation_dyn;
    test<Triangulation_dyn>(D, "inexact dynamic", N);

    typedef CGAL::Epeck_d<CGAL::Dimension_tag<D> > KE;
    typedef CGAL::Delaunay_triangulation<KE> TriangulationE;
    test<TriangulationE>(D, "exact static", N);

    typedef CGAL::Epeck_d<CGAL::Dynamic_dimension_tag> KE_dyn;
    typedef CGAL::Delaunay_triangulation<KE_dyn> TriangulationE_dyn;
    test<TriangulationE_dyn>(D, "exact dynamic", N);
}

int main(int argc, char **argv)
{
    srand(static_cast<unsigned int>(time(nullptr)));
    int N = 10;
    if( argc > 1 )
        N = atoi(argv[1]);
    //go<5>(N);
    go<4>(N);
    go<3>(N);
    go<2>(N);
    go<1>(N);

    cerr << endl;
    return 0;
}

#endif
