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
  ps.add_facet({"a","b","-a","-b","c","d","-c","-d"});
  ps.end_surface();

  CGAL::Path_on_surface<PS> p1(ps); p1.push_back_by_label("a");
  CGAL::Path_on_surface<PS> p2(ps); p2.push_back_by_label({"b","c","a","-c","-b"});
    
  CGAL::Homotopy_tester<PS> ht(ps);

  bool res1=ht.are_freely_homotopic(p1, p2);
  std::cout<<"Paths p1 and p2 "<<(res1?"ARE":"ARE NOT")
           <<" freely homotopic."<<std::endl;

  bool res2=ht.are_base_point_homotopic(p1, p2);
  std::cout<<"Paths p1 and p2 "<<(res2?"ARE":"ARE NOT")
           <<" base point homotopic."<<std::endl;

  return EXIT_SUCCESS;
}

