#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>
#include <cassert>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_descriptor  Dart_descriptor;

int main()
{
  GMap_3 gm;

  // Create one combinatorial hexahedron
  Dart_descriptor d1 = gm.make_combinatorial_hexahedron();

  // Create one square face
  Dart_descriptor d2=gm.make_combinatorial_polygon(4);

  assert(gm.is_insertable_cell_1_between_two_cells_2(d1,d2));

  // Insert the square face as a hole of the face of the hexahedron containing d1
  gm.insert_cell_1_between_two_cells_2(d1, d2);

  // Display the combinatorial map characteristics.
  gm.display_characteristics(std::cout)<<", valid="
                                       <<gm.is_valid()<<std::endl;

  std::size_t nb=0;
  for(Dart_descriptor dh=gm.darts().begin(); dh!=gm.darts().end(); ++dh)
  { if (gm.is_free<2>(dh)) ++nb; }
  std::cout<<"Number of 2-free darts: "<<nb<<std::endl;

  return EXIT_SUCCESS;
}
