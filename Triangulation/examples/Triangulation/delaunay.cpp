#include <CGAL/Cartesian_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <iterator>
#include <iostream>
#include <vector>


typedef CGAL::Cartesian_d<double>           K;
typedef CGAL::Filtered_kernel_d<K>          FK;


typedef CGAL::Triangulation_ds_vertex< void >                   TDS_vertex;
typedef CGAL::Triangulation_vertex< FK, int, TDS_vertex >       Vertex;
typedef CGAL::Triangulation_ds_full_cell
        < void, CGAL::TDS_full_cell_default_storage_policy >    TDS_cell;
typedef CGAL::Triangulation_full_cell< FK, int, TDS_cell >      Cell;
typedef CGAL::Triangulation_data_structure< 
  CGAL::Ambient_dimension< FK::Point_d >::type , Vertex, Cell > TDS;
typedef CGAL::Delaunay_triangulation<FK, TDS>                   T;
//typedef CGAL::Delaunay_triangulation<FK>    T;



int main(int argc, char **argv)
{
    int D = 5; if( argc > 1 )D = atoi(argv[1]);   // space dimension
    int N = 100; if( argc > 2 )N = atoi(argv[2]); // number of points
    CGAL::Timer cost;  // timer

    // Instanciate a random point generator
    CGAL::Random rng;
    typedef CGAL::Random_points_in_cube_d<T::Point> Random_points_iterator;
    Random_points_iterator rand_it(D, 1.0, rng);

    // Generate N random points
    std::vector<T::Point> points;
    CGAL::copy_n(rand_it, N, std::back_inserter(points));

    T t(D);
    assert(t.empty());

    // insert the points in the triangulation
    cost.reset();cost.start();
    std::cout << "  Delaunay triangulation of "<<N<<" points in dim "<<D<< std::flush;
    t.insert(points.begin(), points.end());
    std::cout << " done in "<<cost.time()<<" seconds." << std::endl;
    assert( t.is_valid() );


    return 0;
}
