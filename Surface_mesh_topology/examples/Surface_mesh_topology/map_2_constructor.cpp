#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <iostream>
#include <cstdlib>

using namespace CGAL::Surface_mesh_topology;

typedef Polygonal_schema_with_combinatorial_map<> CMap;
typedef Polygonal_schema_with_generalized_map<>   GMap;

template<typename Map>
void construct_map_from_edges()
{
  Map cm;

  cm.add_facet("a b c -a -b");
  cm.add_facet("-c d e -d -e");

  Path_on_surface<Map> p(cm); p.push_back_by_label("a b d e");

  std::cout<<"Map valid="<<cm.is_valid()<<std::flush;
  cm.display_characteristics(std::cout);
  std::cout<<"; path lenght="<<p.length()
     <<", isclosed? "<<(p.is_closed()?"true":"false")<<std::endl;
}

int main()
{
  construct_map_from_edges<CMap>();
  construct_map_from_edges<GMap>();

  return EXIT_SUCCESS;
}

