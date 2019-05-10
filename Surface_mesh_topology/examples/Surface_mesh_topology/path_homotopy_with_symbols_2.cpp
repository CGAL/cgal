#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_2_incremental_builder.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Homotopy_tester.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<2> CMap;

int main()
{
  CMap cm;
  CGAL::Combinatorial_map_2_incremental_builder<CMap> b(cm);

  b.begin_surface();

  b.add_facet("a b -a c"); // First facet, giving directly its sequence of edges
  b.add_facet("d -c e -b"); // Second facet

  b.begin_facet(); // Third facet
  b.add_edges_to_facet("f"); // Here, each edge is added one at a time
  b.add_edges_to_facet("-d");
  b.add_edges_to_facet("-f");
  b.add_edges_to_facet("-e");
  b.end_facet();
    
  b.end_surface();

  std::cout<<"Number of cells of the combinatorial maps: ";
  cm.display_characteristics(std::cout)<<std::endl;
  
  CGAL::Path_on_surface<CMap> p=b.create_path("a b -a e -b d");  
  CGAL::Homotopy_tester<CMap> smct(cm);

  bool res=smct.is_contractible(p);
  std::cout<<"Path "<<(res?"IS":"IS NOT")<<" contractible."<<std::endl;

  return EXIT_SUCCESS;
}

