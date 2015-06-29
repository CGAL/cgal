#include <fstream>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <boost/foreach.hpp>

struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  // The point type that will be injected here will be
  // CGAL::Exact_predicates_inexact_constructions_kernel::Point_3
  template <typename Point>
  bool operator()(const Point& p, const Point& q, const Point& r) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return false;
    }
    double d  = sqrt(squared_distance(p,q));
    if(d>bound) return true;
    d += sqrt(squared_distance(p,r)) ;
    if(d>bound) return true;
    d+= sqrt(squared_distance(q,r));
    return d>bound;
  }
};


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Advancing_front_surface_reconstruction<CGAL::Default, Perimeter> Reconstruction;
typedef Reconstruction::Triangulation_3 Triangulation_3;
typedef Reconstruction::Outlier_range Outlier_range;
typedef Reconstruction::Boundary_range Boundary_range;
typedef Reconstruction::Vertex_on_boundary_range Vertex_on_boundary_range;
typedef Reconstruction::Vertex_handle Vertex_handle;
typedef K::Point_3 Point_3;


int main(int argc, char* argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/half.xyz");
  std::istream_iterator<Point_3> begin(in);
  std::istream_iterator<Point_3> end;

  Perimeter perimeter(0.5);
  Triangulation_3 dt(begin, end);
  Reconstruction reconstruction(dt, perimeter);

  reconstruction.run();

  std::cout << reconstruction.number_of_outliers() << " outliers:\n" << std::endl;
  BOOST_FOREACH(const Point_3& p, reconstruction.outliers()){
    std::cout << p << std::endl;
  }

  std::cout << "Boundaries:" << std::endl ;
  BOOST_FOREACH(const Vertex_on_boundary_range & vobr, reconstruction.boundaries()){
    std::cout << "boundary\n";
    // As we use BOOST_FOREACH we do not use the type Boundary_range
    BOOST_FOREACH(Vertex_handle v, vobr){
      std::cout << v->point() << std::endl;
    }
  }

  return 0;
}
