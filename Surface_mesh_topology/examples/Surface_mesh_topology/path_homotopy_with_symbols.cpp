#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_2_incremental_builder.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_curve_topology.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<2> CMap;

int main()
{
  CMap cm;
  CGAL::Combinatorial_map_2_incremental_builder<CMap> b(cm);

  b.begin_surface();
  b.add_facet("a b -a -b c d -c -d");
  b.end_surface();

  CGAL::Path_on_surface<CMap> p1=b.create_path("a");
  CGAL::Path_on_surface<CMap> p2=b.create_path("b c a -c -b");
    
  CGAL::Surface_mesh_curve_topology<CMap> smct(cm);

  bool res1=smct.are_freely_homotopic(p1, p2);
  std::cout<<"Paths p1 and p2 "<<(res1?"ARE":"ARE NOT")
           <<" freely homotopic."<<std::endl;

  bool res2=smct.are_base_point_homotopic(p1, p2);
  std::cout<<"Paths p1 and p2 "<<(res2?"ARE":"ARE NOT")
           <<" base point homotopic."<<std::endl;

  return EXIT_SUCCESS;
}

