#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<2> GMap_2;
typedef GMap_2::Dart_descriptor Dart_descriptor;

int main()
{
  GMap_2 gm;
  Dart_descriptor d=gm.make_combinatorial_polygon(4);
  gm.sew<2>(d, gm.alpha<1,0,1,0>(d));

  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
