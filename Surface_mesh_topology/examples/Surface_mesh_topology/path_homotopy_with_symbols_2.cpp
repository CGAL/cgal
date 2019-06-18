#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Homotopy_tester.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Polygonal_schema_with_combinatorial_map<> PS;

int main()
{
  PS ps;

  ps.begin_surface();

  ps.add_facet("a b -a c"); // First facet, giving directly its sequence of edges
  ps.add_facet("d -c e -b"); // Second facet

  ps.begin_facet(); // Third facet
  ps.add_edges_to_facet("f"); // Here, each edge is added one at a time
  ps.add_edges_to_facet("-d");
  ps.add_edges_to_facet("-f");
  ps.add_edges_to_facet("-e");
  ps.end_facet();
    
  ps.end_surface();

  std::cout<<"Number of cells of the combinatorial maps: ";
  ps.display_characteristics(std::cout)<<std::endl;
  
  CGAL::Path_on_surface<PS> p(ps);
  p.push_back_by_label("a b -a e -b d");  

  CGAL::Homotopy_tester<PS> ht(ps);
  bool res=ht.is_contractible(p);
  std::cout<<"Path "<<(res?"IS":"IS NOT")<<" contractible."<<std::endl;

  return EXIT_SUCCESS;
}

