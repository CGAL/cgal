#include <CGAL/Epick_d.h>
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

    typedef CGAL::Random_points_in_cube_d<Point> Random_points_iterator;
    CGAL::Timer cost;  // timer

    DT dt(d);
    assert(dt.empty());

    std::vector<Point> points;
    CGAL::Random rng;
    Random_points_iterator rand_it(d, 2.0, rng);
    CGAL::cpp11::copy_n(rand_it, N, std::back_inserter(points));
    cost.reset();cost.start();
    std::cout << "  Delaunay triangulation of "<<N<<" points in dim "<<d<< std::flush;
    dt.insert(points.begin(), points.end());
    std::cout << " done in "<<cost.time()<<" seconds." << std::endl;
    std::size_t nbfc= dt.number_of_finite_full_cells();
    std::size_t nbc= dt.number_of_full_cells();
    std::cout << dt.number_of_vertices() << " vertices, " 
              << nbfc << " finite simplices and " 
              << (nbc-nbfc) << " convex hull Facets."
              << std::endl;
}

template< int D >
void go(const int N)
{
    typedef CGAL::Epick_d<CGAL::Dimension_tag<D> > K;
    //typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag> K;
    typedef CGAL::Delaunay_triangulation<K> Triangulation;
    test<Triangulation>(D, "static", N);
}

int main(int argc, char **argv)
{
    srand(static_cast<unsigned int>(time(NULL)));
    int N = 100; if( argc > 1 ) N = atoi(argv[1]);
    go<2>(N);
    go<3>(N);
    go<4>(N);
    go<5>(N);
    go<6>(N);
    go<7>(N);
    go<8>(N);


    return 0;
}
