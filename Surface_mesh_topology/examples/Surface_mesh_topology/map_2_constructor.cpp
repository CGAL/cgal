#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_incremental_builder.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<2> CMap;

int main()
{
  CMap cm;
  CGAL::Combinatorial_map_incremental_builder<CMap> b(cm);

  b.begin_surface(0,0,0);
  b.add_facet("a b c -a -b");
  b.add_facet("-c d e -d -e");
  b.end_surface();
  
  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;

  Path p=b.create_path("a b -a");
  
  
  return EXIT_SUCCESS;
}

