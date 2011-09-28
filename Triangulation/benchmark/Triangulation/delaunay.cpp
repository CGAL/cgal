#include <CGAL/Cartesian_d.h>
//#include <CGAL/Simple_cartesian_d.h>
//#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Timer.h>
#include <CGAL/algorithm.h>
#include <vector>
#include <string>
#include <fstream>
#include <cstdlib>
#include <algorithm>


template<typename DT>
void test(const int d, const std::string & type, const int N)
{
    typedef typename DT::Vertex Vertex;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Full_cell Full_cell;
    typedef typename DT::Full_cell_handle Full_cell_handle;
    typedef typename DT::Facet Facet;
    typedef typename DT::Point Point;
    typedef typename DT::Geom_traits::RT RT;
    typedef typename DT::Finite_full_cell_const_iterator Finite_full_cell_const_iterator;

    typedef CGAL::Random_points_in_iso_box_d<Point> Random_points_iterator;
    CGAL::Timer cost;  // timer

    DT dt(d);
    assert(dt.empty());

    std::vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 2.0, rng);
    CGAL::copy_n(rand_it, N, std::back_inserter(points));
    cost.reset();cost.start();
    std::cout << "  Delaunay triangulation of "<<N<<" points in dim "<<d<< std::flush;
    dt.insert(points.begin(), points.end());
    std::cout << " done in "<<cost.time()<<" seconds." << std::endl;
    int nbfc= dt.number_of_finite_full_cells();
    int nbc= dt.number_of_full_cells();
    std::cout << dt.number_of_vertices() << " vertices, " 
	      << nbfc << " finite simplices and " 
	      << (nbc-nbfc) << " convex hull Facets."
	      << std::endl;
}

template< int D, typename RT >
void go(const int N)
{
    //typedef CGAL::Simple_cartesian_d<RT, D> K;
    typedef CGAL::Cartesian_d<RT> K;
    //typedef CGAL::Filtered_kernel_d<K> FK;
    typedef CGAL::Delaunay_triangulation<K> Triangulation;
    test<Triangulation>(D, "static", N);
}

int main(int argc, char **argv)
{
    srand48(time(NULL));
    int N = 100; if( argc > 1 ) N = atoi(argv[1]);
    go<2, double>(N);
    go<3, double>(N);
    go<4, double>(N);
    go<5, double>(N);
    go<6, double>(N);
    go<7, double>(N);


    return 0;
}
