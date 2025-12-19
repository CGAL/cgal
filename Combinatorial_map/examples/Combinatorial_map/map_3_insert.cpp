#include <CGAL/Combinatorial_map.h>
#include <iostream>
#include <cstdlib>
#include <cassert>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_descriptor    Dart_descriptor;

int main()
{
  CMap_3 cm;

  // Create one combinatorial hexahedron
  Dart_descriptor d1 = cm.make_combinatorial_hexahedron();

  // Create one square face
  Dart_descriptor d2=cm.make_combinatorial_polygon(4);

  assert(cm.is_insertable_cell_1_between_two_cells_2(d1,d2));

  // Insert the square face as a hole of the face of the hexahedron containing d1
  cm.insert_cell_1_between_two_cells_2(d1, d2);

  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout)<<", valid="
                                       <<cm.is_valid()<<std::endl;

  std::size_t nb=0;
  for(Dart_descriptor dh=cm.darts().begin(); dh!=cm.darts().end(); ++dh)
  { if (cm.is_free<2>(dh)) ++nb; }
  std::cout<<"Number of 2-free darts: "<<nb<<std::endl;

  return EXIT_SUCCESS;
}
