#include <fstream>
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>

struct Perimeter {

  double bound;

  Perimeter(double bound)
    : bound(bound)
  {}

  template <typename AdvancingFront, typename Cell_handle>
  double operator() (const AdvancingFront& adv, Cell_handle& c,
                     const int& index) const
  {
    // bound == 0 is better than bound < infinity
    // as it avoids the distance computations
    if(bound == 0){
      return adv.smallest_radius_delaunay_sphere (c, index);
    }

    // If perimeter > bound, return infinity so that facet is not used
    double d  = 0;
    d = sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                              c->vertex((index+2)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+2)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();
    d += sqrt(squared_distance(c->vertex((index+1)%4)->point(),
                               c->vertex((index+3)%4)->point()));
    if(d>bound) return adv.infinity();

    // Otherwise, return usual priority value: smallest radius of
    // delaunay sphere
    return adv.smallest_radius_delaunay_sphere (c, index);
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
  for(const Point_3& p : reconstruction.outliers()){
    std::cout << p << std::endl;
  }

  std::cout << "Boundaries:" << std::endl ;
  for(const Vertex_on_boundary_range& vobr : reconstruction.boundaries()){
    std::cout << "boundary\n";
    // As we use range-base loop we do not use the type Boundary_range
    for(Vertex_handle v : vobr){
      std::cout << v->point() << std::endl;
    }
  }

  return 0;
}
