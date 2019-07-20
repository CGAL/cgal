#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<4> GMap_4;
typedef GMap_4::Dart_handle Dart_handle;

int main()
{
  GMap_4 gm;
  Dart_handle d1 = gm.make_combinatorial_tetrahedron();
  Dart_handle d2 = gm.make_combinatorial_tetrahedron();

  CGAL_assertion(gm.is_valid());

  gm.sew<4>(d1,d2);

  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
