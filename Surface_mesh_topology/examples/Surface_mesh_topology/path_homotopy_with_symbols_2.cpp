#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Surface_mesh_topology::Polygonal_schema_with_combinatorial_map<> PS;

int main()
{
  PS ps;

  ps.add_facet("a b -a c"); // First facet, giving directly its sequence of edges
  ps.add_facet("d -c e -b"); // Second facet

  ps.init_facet(); // Third facet
  ps.add_edges_to_facet("f"); // Here, each edge is added one at a time
  ps.add_edges_to_facet("-d");
  ps.add_edges_to_facet("-f");
  ps.add_edges_to_facet("-e");
  ps.finish_facet();

  ps.perforate_facet("f");

  std::cout<<"Number of cells of the combinatorial maps: ";
  ps.display_characteristics(std::cout)<<std::endl;

  CGAL::Surface_mesh_topology::Path_on_surface<PS> p(ps);
  p.push_back_by_label("a b -a e -b d");

  CGAL::Surface_mesh_topology::Curves_on_surface_topology<PS> cst(ps);
  bool res=cst.is_contractible(p);
  std::cout<<"Path "<<(res?"IS":"IS NOT")<<" contractible."<<std::endl;

  return EXIT_SUCCESS;
}

