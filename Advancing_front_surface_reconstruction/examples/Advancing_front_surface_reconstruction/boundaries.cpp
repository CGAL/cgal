#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Advancing_front_surface_reconstruction<K> Reconstruction;
typedef Reconstruction::Outlier_range Outlier_range;
typedef Reconstruction::Boundary_range Boundary_range;
typedef Reconstruction::Vertex_on_boundary_range Vertex_on_boundary_range;
typedef Reconstruction::Vertex_handle Vertex_handle;
typedef Kernel::Point_3 Point_3;

int main()
{
  std::istream_iterator<Point_3> begin(std::cin);
  std::istream_iterator<Point_3> end;

  Triangulation_3 dt(begin,end); 
  Reconstruction reconstruction(dt);

  reconstruction.run();
                
  std::cout << "Outliers:\n";
  Outlier_range::iterator b,e;
  // instead of boost::tie we might use BOOST_FOREACH
  for(boost::tie(b,e) = reconstruction.outliers()){
    Point_3 p = *b;
    std::cout << p << std::endl;
  }
  
  std::cout << "Boundaries:\n";
  BOOST_FOREACH(Vertex_on_boundary_range it, reconstruction.boundaries()){
    std::cout << "boundary\n";
    // As we use BOOST_FOREACH we do not use the type Boundary_range
    BOOST_FOREACH(Vertex_handle v, *it){
      std::cout << v->point() << std::endl;
    }
  }  
  return 0;
}
