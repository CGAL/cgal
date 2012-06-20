#include <CGAL/Cartesian_d.h>
//#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/point_generators_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/algorithm.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <iterator>
#include <iostream>
#include <vector>
const int D=5;

typedef CGAL::Cartesian_d<double>                              K;//D;
//typedef CGAL::Filtered_kernel_d<KD>                            K;
typedef CGAL::Triangulation_ds_vertex< void >                  TDS_vertex;
typedef CGAL::Triangulation_vertex< K, int, TDS_vertex >       Vertex;
typedef CGAL::Triangulation_ds_full_cell
        < void, CGAL::TDS_full_cell_default_storage_policy >    TDS_cell;
typedef CGAL::Triangulation_full_cell< K, int, TDS_cell >       Cell;
typedef CGAL::Triangulation_data_structure< 
   CGAL::Dimension_tag<D> , Vertex, Cell >  TDS;
typedef CGAL::Delaunay_triangulation<K,TDS>                    T;


int main(int argc, char **argv)
{
  int N = 100; if( argc > 2 )N = atoi(argv[1]); // number of points
  CGAL::Timer cost;  // timer
  
  // Instanciate a random point generator
  CGAL::Random rng(0);
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

  // insert with special operations in conflict zone and new created cells
  cost.reset();
  std::cout << "  adding "<<N<<" other points "<< std::endl;
  for(int i=0; i<N; ++i){
    T::Vertex_handle v;
    T::Face	f(t.ambient_dimension()); 
    T::Facet	ft; 
    T::Full_cell_handle c; 
    T::Locate_type	lt;
    typedef std::vector<T::Full_cell_handle> Full_cells; 
    Full_cells zone, new_full_cells; 
    std::back_insert_iterator<Full_cells> out(zone); 
    c = t.locate(*++rand_it, lt, f, ft, v);
    // previously inserted vertex v is used as hint for point location (if defined)
    T::Facet ftc = t.compute_conflict_zone(*rand_it, c, out); 
    std::cout<<i<<"     conflict zone of size "<<zone.size()<<" -> "<<std::flush;
    out = std::back_inserter(new_full_cells);
    std::cout<<" locate type "<<lt<<"; "<<std::flush;
    v = t.insert_in_hole(*rand_it, zone.begin(), zone.end(), ftc, out);
    std::cout<<new_full_cells.size()<<" new cells"<<std::endl;
    // for (Full_cells::iterator it=new_full_cells.begin();
    //	 it!=new_full_cells.end(); ++it) (*it)->data() = zone.size();
  }
  std::cout << " done in "<<cost.time()<<" seconds." << std::endl;



  return 0;
}
