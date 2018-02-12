#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

typedef CGAL::Surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> SM;

int main()
{
  SM sm;
  std::ifstream is("data/self-intersection.off");
  is>>sm;
  
  bool is_ok =
  CGAL::Polygon_mesh_processing::remove_self_intersections(sm, 
                                                           CGAL::Polygon_mesh_processing::parameters::all_default(), 
                                                           7,
                                                           true);
  assert(is_ok);
  return 0;
}