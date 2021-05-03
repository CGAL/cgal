#include <CGAL/config.h>
#if defined(BOOST_GCC) && (__GNUC__ <= 4) && (__GNUC_MINOR__ < 4)
#include <iostream>
int main()
{
  std::cerr << "NOTICE: This test requires G++ >= 4.4, and will not be compiled." << std::endl;
}

#else
#include <CGAL/Epick_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Triangulation.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <CGAL/IO/Triangulation_off_ostream.h>

using namespace std;

template<typename T>
void test(const int d, const string & type, int N)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that T contains sub-types, before
    // instanciating the function.
    typedef typename T::Full_cell_handle Full_cell_handle;
    typedef typename T::Point Point;
    typedef typename T::Geom_traits::RT RT;
    typedef typename T::Finite_vertex_const_iterator Finite_vertex_const_iterator;
    typedef typename T::Finite_facet_iterator Finite_facet_iterator;
    typedef typename T::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;

    typedef CGAL::Random_points_in_cube_d<Point> Random_points_iterator;

    T tri(d);
    cerr << "\nChecking Triangulation of (" << type << d << ") dimension "
        << tri.maximal_dimension();
    assert(tri.empty());

    vector<RT> coords(d);
    vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 1.0, rng);
    std::copy_n(rand_it, N, std::back_inserter(points));

    cerr << '\n' << points.size() << " points in the grid.";

    tri.insert(points.begin(), points.end());
    assert( tri.is_valid() );

    cerr << "\nTraversing finite full_cells... ";
    size_t nbfs(0), nbis(0);
    Finite_full_cell_const_iterator fsit = tri.finite_full_cells_begin();
    while( fsit != tri.finite_full_cells_end() )
    {
        ++fsit, ++nbfs;
    }
    cerr << nbfs << " + ";
    vector<Full_cell_handle> infinite_full_cells;
    tri.tds().incident_full_cells(tri.infinite_vertex(), std::back_inserter(infinite_full_cells));
    nbis = infinite_full_cells.size();
    cerr << nbis << " = " << (nbis+nbfs)
    << " = " << tri.number_of_full_cells();
    assert(nbfs + nbis == tri.number_of_full_cells());

    cerr << "\nTraversing finite facets... ";
    size_t nbff(0);
    Finite_facet_iterator ffit = tri.finite_facets_begin();
    while( ffit != tri.finite_facets_end() )
    {
        ++ffit, ++nbff;
    }
    cerr << nbff << " finite facets";

    cerr << "\nTraversing finite vertices... ";
    size_t nbfv(0);
    Finite_vertex_const_iterator fvit = tri.finite_vertices_begin();
    while( fvit != tri.finite_vertices_end() )
    {
        ++fvit, ++nbfv;
    }
    cerr << nbfv << " finite vertices (should be " << tri.number_of_vertices() << ").";
    assert(nbfv == tri.number_of_vertices());

    // TEST Copy Constructor
    T tri2(tri);
    assert( tri2.is_valid() );
    assert( tri.current_dimension() == tri2.current_dimension() );
    assert( tri.maximal_dimension() == tri2.maximal_dimension() );
    assert( tri.number_of_vertices() == tri2.number_of_vertices() );
    assert( tri.number_of_full_cells() == tri2.number_of_full_cells() );

    std::stringstream buffer;
    buffer << tri;

    // CLEAR
    tri.clear();
    assert(-1==tri.current_dimension());
    assert(tri.empty());
    assert( tri.is_valid() );

    buffer >> tri;
    assert( tri.current_dimension() == tri2.current_dimension() );
    assert( tri.maximal_dimension() == tri2.maximal_dimension() );
    assert( tri.number_of_vertices() == tri2.number_of_vertices() );
    assert( tri.number_of_full_cells() == tri2.number_of_full_cells() );

    std::ofstream ofs("tri", std::ios::binary);
    ofs << tri;
    ofs.close();

    std::ifstream ifs("tri", std::ios::binary);
    ifs >> tri2;
    ifs.close();
    assert( tri.current_dimension() == tri2.current_dimension() );
    assert( tri.maximal_dimension() == tri2.maximal_dimension() );
    assert( tri.number_of_vertices() == tri2.number_of_vertices() );
    assert( tri.number_of_full_cells() == tri2.number_of_full_cells() );
}

/*#define test_static(DIM) {  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<DIM>> T; \
        test<T>(DIM, string("static")+string(#DIM)); }
#define test_dyn(DIM) { \
    typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag> T; \
        test<T>(DIM, string("dynamic")+string(#DIM)) ;}

#define test_mirror_static(DIM) {  \
    typedef CGAL::Triangulation_ds_full_cell<void, CGAL::T_full_cell_mirror_storage_policy> My_ds_full_cell;  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dimension_tag<DIM>, \
                                              CGAL::Triangulation_ds_vertex<>, \
                                              My_ds_full_cell> My_tds;  \
        test<My_tds>(DIM, string("mirror&static")+string(#DIM)); }
#define test_mirror_dyn(DIM) { \
    typedef CGAL::Triangulation_ds_full_cell<void, CGAL::T_full_cell_mirror_storage_policy> My_ds_full_cell;  \
    typedef CGAL::Triangulation_data_structure<CGAL::Dynamic_dimension_tag, \
                                              CGAL::Triangulation_ds_vertex<>, \
                                              My_ds_full_cell> My_tds;  \
        test<My_tds>(DIM, string("mirror&dynamic")+string(#DIM)) ;}
*/

template< int D >
void go(int N)
{
    typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > K;
    //typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
    typedef CGAL::Triangulation<K> Triangulation;
    test<Triangulation>(D, "static", N);
    //test<Triangulation>(D, "dynamic", N);
}

int main(int argc, char **argv)
{
    srand(static_cast<unsigned int>(time(nullptr)));
    int N = 1000;
    if( argc > 1 )
        N = atoi(argv[1]);
    //go<5>(N);
    go<4>(N);
    go<3>(N);
    go<2>(N);
    go<1>(N);

    cerr << std::endl;
    return 0;
}

#endif
