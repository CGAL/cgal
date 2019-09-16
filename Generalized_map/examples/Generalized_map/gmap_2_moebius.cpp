#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<2> GMap_2;
typedef GMap_2::Dart_handle Dart_handle;

int main()
{
  GMap_2 gm;
  Dart_handle dh=gm.make_combinatorial_polygon(4);
  gm.sew<2>(dh, gm.alpha<1,0,1,0>(dh));

  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
