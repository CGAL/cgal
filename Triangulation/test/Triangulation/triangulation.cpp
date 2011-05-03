#include <CGAL/point_generators_d.h>
//#include <CGAL/Simple_cartesian_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Pure_complex.h>
#include <CGAL/algorithm.h>
#include <tilted_grid.h>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

template<typename PC>
void test(const int d, const string & type, int N)
{
    // we must write 'typename' below, because we are in a template-function,
    // so the parser has no way to know that PC contains sub-types, before
    // instanciating the function.
    typedef typename PC::Vertex Vertex;
    typedef typename PC::Vertex_handle Vertex_handle;
    typedef typename PC::Simplex Simplex;
    typedef typename PC::Simplex_handle Simplex_handle;
    typedef typename PC::Facet Facet;
    typedef typename PC::Point Point;
    typedef typename PC::Geom_traits::RT RT;
    typedef typename PC::Finite_simplex_const_iterator Finite_simplex_const_iterator;

    typedef CGAL::Random_points_in_iso_box_d<Point> Random_points_iterator;

    PC pc(d);
    cerr << "\nChecking Pure_complex of (" << type << d << ") dimension "
        << pc.ambient_dimension();
    assert(pc.empty());

    vector<RT> coords(d);
	vector<Point> points;
    CGAL::Random rng;
	Random_points_iterator rand_it(d, 1.0, rng);
	CGAL::copy_n(rand_it, N, std::back_inserter(points));

    cerr << '\n' << points.size() << " points in the grid.";

    pc.insert(points.begin(), points.end());
    assert( pc.is_valid() );

    cerr << "\nTraversing finite simplices... ";
    size_t nbfs(0), nbis(0);
    Finite_simplex_const_iterator fsit = pc.finite_simplices_begin();
    while( fsit != pc.finite_simplices_end() )
    {
        Point c = fsit->circumcenter();
        ++fsit, ++nbfs;
    }
    cerr << nbfs << " + ";
    vector<Simplex_handle> infinite_simplices;
    pc.gather_incident_simplices(pc.infinite_vertex(), std::back_inserter(infinite_simplices));
    nbis = infinite_simplices.size();
    cerr << nbis << " = " << (nbis+nbfs)
    << " = " << pc.number_of_simplices();

    // CLEAR
    pc.clear();
    assert(-1==pc.current_dimension());
    assert(pc.empty());
    assert( pc.is_valid() );
}

/*#define test_static(DIM) {  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dimension_tag<DIM>> PC; \
        test<PC>(DIM, string("static")+string(#DIM)); }
#define test_dyn(DIM) { \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dynamic_dimension_tag> PC; \
        test<PC>(DIM, string("dynamic")+string(#DIM)) ;}

#define test_mirror_static(DIM) {  \
    typedef CGAL::Pure_complex_ds_simplex<void, CGAL::PC_simplex_mirror_storage_policy> My_ds_simplex;  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dimension_tag<DIM>, \
                                              CGAL::Pure_complex_ds_vertex<>, \
                                              My_ds_simplex> My_pcds;  \
        test<My_pcds>(DIM, string("mirror&static")+string(#DIM)); }
#define test_mirror_dyn(DIM) { \
    typedef CGAL::Pure_complex_ds_simplex<void, CGAL::PC_simplex_mirror_storage_policy> My_ds_simplex;  \
    typedef CGAL::Pure_complex_data_structure<CGAL::Dynamic_dimension_tag, \
                                              CGAL::Pure_complex_ds_vertex<>, \
                                              My_ds_simplex> My_pcds;  \
        test<My_pcds>(DIM, string("mirror&dynamic")+string(#DIM)) ;}
*/

template< int D >
void go(int N)
{
    typedef double RT;
    //typedef CGAL::Simple_cartesian_d<RT, D> K;
    typedef CGAL::Cartesian_d<RT> K;
    typedef CGAL::Filtered_kernel_d<K> FK;
    typedef CGAL::Pure_complex<FK> Triangulation;
    //test<Triangulation>(D, "static", N);
    test<Triangulation>(D, "dynamic", N);
}

int main(int argc, char **argv)
{
    srand48(time(NULL));
    int N = 1000;
    if( argc > 1 )
        N = atoi(argv[1]);
    // go<5>(N);
     go<3>(N);
    // go<2>(N);
    // go<1>(N);

    cerr << std::endl;
    return 0;
}
