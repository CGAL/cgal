#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/repair.h>

typedef CGAL::Surface_mesh<CGAL::Exact_predicates_inexact_constructions_kernel::Point_3> SM;

int main()
{
  namespace PMP = CGAL::Polygon_mesh_processing;
  SM sm;
  std::ifstream is("data/self_intersection.off");
  is>>sm;
  
  bool is_ok =
  PMP::remove_self_intersections(sm, 
                                                           PMP::parameters::all_default(), 
                                                           7,
                                                           false);
  assert(is_ok);
  std::ofstream os("/home/gimeno/Data/tmp/debug.off");
     os << sm;
  std::vector< std::pair<boost::graph_traits<SM>::face_descriptor,
      boost::graph_traits<SM>::face_descriptor> > sis;
  PMP::self_intersections(sm,
                          std::back_inserter(sis),
                          PMP::parameters::all_default());
  assert(sis.empty());
  return 0;
}