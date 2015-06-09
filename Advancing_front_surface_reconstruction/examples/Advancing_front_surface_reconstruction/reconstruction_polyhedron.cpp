

#include <iostream>
#include <algorithm>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Advancing_front_surface_reconstruction.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3  Point_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;


int main()
{
  std::vector<Point_3> points;
  Polyhedron P;

  std::copy(std::istream_iterator<Point_3>(std::cin), 
            std::istream_iterator<Point_3>(), 
            std::back_inserter(points));
  
  CGAL::advancing_front_surface_reconstructionP(points.begin(),
                                               points.end(),
                                               P);

  std::cout << P  << std::endl;

  return 0;
}
