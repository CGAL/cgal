#include <CGAL/Combinatorial_map.h>
#include <CGAL/Generalized_map.h>
#include <CGAL/Combinatorial_map_2_incremental_builder.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<2> CMap;
typedef CGAL::Generalized_map<2>   GMap;

template<typename Map>
void construct_map_from_edges()
{
  Map cm;
  CGAL::Combinatorial_map_2_incremental_builder<Map> b(cm);

  b.begin_surface();
  b.add_facet("a b c -a -b");
  b.add_facet("-c d e -d -e");
  b.end_surface();

  CGAL::Path_on_surface<Map> p=b.create_path("a b d e");

  std::cout<<"Map valid="<<cm.is_valid()<<std::flush;
  cm.display_characteristics(std::cout);
  std::cout<<"; path lenght="<<p.length()
           <<", isclosed? "<<(p.is_closed()?"true":"false")<<std::endl;
}

template<typename Map>
void construct_map_from_vertices()
{
  Map cm;
  CGAL::Combinatorial_map_2_incremental_builder<Map> b(cm);

  b.begin_surface();

  b.begin_facet();
  b.add_vertex_to_facet(0); b.add_vertex_to_facet(1); b.add_vertex_to_facet(2);
  b.end_facet();

  b.begin_facet();
  b.add_vertex_to_facet(1); b.add_vertex_to_facet(0); b.add_vertex_to_facet(3);
  b.end_facet();

  b.begin_facet();
  b.add_vertex_to_facet(2); b.add_vertex_to_facet(3); b.add_vertex_to_facet(0);
  b.end_facet();

  b.begin_facet();
  b.add_vertex_to_facet(3); b.add_vertex_to_facet(2); b.add_vertex_to_facet(1);
  b.end_facet();

  b.begin_path();
  b.add_vertex_to_path(0);
  b.add_vertex_to_path(3);
  b.add_vertex_to_path(2);
  b.add_vertex_to_path(1);
  b.add_vertex_to_path(0);
  CGAL::Path_on_surface<Map> p=b.end_path();

  std::cout<<"Map valid="<<cm.is_valid()<<std::flush;
  cm.display_characteristics(std::cout);
  std::cout<<"; path lenght="<<p.length()
           <<", isclosed? "<<(p.is_closed()?"true":"false")<<std::endl;
}

int main()
{
  construct_map_from_vertices<CMap>();
  construct_map_from_edges<CMap>();
  construct_map_from_vertices<GMap>();
  construct_map_from_edges<GMap>();
  return EXIT_SUCCESS;
}

