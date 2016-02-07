#if defined(__GNUC__) && defined(__GNUC_MINOR__) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)

#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else

#include <CGAL/Epick_d.h>
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

    typedef CGAL::Random_points_in_cube_d<Point> Random_points_iterator;

    DC pc(d);
    cerr << "\nBuilding Delaunay triangulation of (" << type << d << ") dimension with " << N << " points";
    assert(pc.empty());

    vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 2.0, rng);
    //CGAL::cpp11::copy_n(rand_it, N, back_inserter(points));
    
    vector<int> coords(d);
    for( int i = 0; i < N; ++i )
    {
        for( int j = 0; j < d; ++j )
            coords[j] = rand() % 100000;
        points.push_back(Point(d, coords.begin(), coords.end()));
    }
    pc.insert(points.begin(),  points.end());
    cerr << "\nChecking topology and geometry...";
    assert( pc.is_valid() );

    cerr << "\nTraversing finite full_cells... ";
    size_t nbfs(0), nbis(0);
    Finite_full_cell_const_iterator fsit = pc.finite_full_cells_begin();
    while( fsit != pc.finite_full_cells_end() )
        ++fsit, ++nbfs;
    cerr << nbfs << " + ";
    vector<Full_cell_handle> infinite_full_cells;
    pc.tds().incident_full_cells(pc.infinite_vertex(), back_inserter(infinite_full_cells));
    nbis = infinite_full_cells.size();
    cerr << nbis << " = " << (nbis+nbfs)
    << " = " << pc.number_of_full_cells();
    cerr << "\nThe triangulation has current dimension " << pc.current_dimension();
    CGAL_assertion( pc.number_of_full_cells() == nbis+nbfs);

    cerr << "\nTraversing finite vertices... ";
    size_t nbfv(0);
    Finite_vertex_iterator fvit = pc.finite_vertices_begin();
    while( fvit != pc.finite_vertices_end() )
        ++fvit, ++nbfv;
    cerr << nbfv <<endl;

    // Count convex hull vertices:
    if( pc.maximal_dimension() > 1 )
    {
        typedef vector<Face> Faces;
        Faces edges;
        back_insert_iterator<Faces> out(edges);
        pc.tds().incident_faces(pc.infinite_vertex(), 1, out);
        cout << "\nThere are " << edges.size() << " vertices on the convex hull.";
        edges.clear();
    }
    else // pc.maximal_dimension() == 1
    {
        typedef vector<Full_cell_handle> Cells;
        Cells cells;
        back_insert_iterator<Cells> out(cells);
        pc.tds().incident_full_cells(pc.infinite_vertex(), out);
        cout << "\nThere are " << cells.size() << " vertices on the convex hull.";
        cells.clear();
    }

    // Remove all !
    cerr << "\nBefore removal: " << pc.number_of_vertices() << " vertices. After: ";
    random_shuffle(points.begin(),  points.end());
    pc.remove(points.begin(),  points.end());
    assert( pc.is_valid() );
    cerr << pc.number_of_vertices() << " vertices.";
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
    //typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> FK;
    typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > FK;
    typedef CGAL::Delaunay_triangulation<FK> Triangulation;
    //test<Triangulation>(D, "dynamic", N);
    test<Triangulation>(D, "static", N);
}

int main(int argc, char **argv)
{
    srand(static_cast<unsigned int>(time(NULL)));
    int N = 100;
    if( argc > 1 )
        N = atoi(argv[1]);
     go<5>(N);
     go<4>(N);
     go<3>(N);
     go<2>(N);
     go<1>(N);

    cerr << endl;
    return 0;
}

#endif
