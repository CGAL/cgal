#include <CGAL/Combinatorial_map.h>
#include <iostream>
#include <cstdlib>
#include <cassert>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_descriptor    Dart_descriptor;

int main()
{
  CMap_3 cm;

  // Create one combinatorial hexahedron.
  Dart_descriptor d1 = cm.make_combinatorial_hexahedron();

  // Add two edges along two opposite facets.
  assert( cm.is_insertable_cell_1_in_cell_2
          (cm.beta(d1,1),cm.beta(d1,0)) );

  cm.insert_cell_1_in_cell_2(cm.beta(d1,1), cm.beta(d1,0));
  assert( cm.is_valid() );

  Dart_descriptor d2=cm.beta(d1,2,1,1,2);

  assert( cm.is_insertable_cell_1_in_cell_2
          (d2,cm.beta(d2,1,1)) );

  cm.insert_cell_1_in_cell_2(d2, cm.beta(d2,1,1));
  assert( cm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges
  // of the hexahedron.
  std::vector<Dart_descriptor> path;
  path.push_back(cm.beta(d1,1));
  path.push_back(cm.beta(d1,0,2,1));
  path.push_back(cm.beta(d2,0));
  path.push_back(cm.beta(d2,2,1));

  assert( (cm.is_insertable_cell_2_in_cell_3
                               (path.begin(),path.end())) );

  Dart_descriptor d3=cm.insert_cell_2_in_cell_3(path.begin(),path.end());
  assert( cm.is_valid() );

  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid=" <<
    cm.is_valid() << std::endl;

  // We use the removal operations to get back to the initial hexahedron.
  assert( (cm.is_removable<2>(d3)) );
  cm.remove_cell<2>(d3);
  assert( cm.is_valid() );

  assert( (cm.is_removable<1>(cm.beta(d1,1))) );
  cm.remove_cell<1>(cm.beta(d1,1));
  assert( cm.is_valid() );

  assert( (cm.is_removable<1>(cm.beta(d2,0))) );
  cm.remove_cell<1>(cm.beta(d2,0));
  assert( cm.is_valid() );

  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid="
                                        << cm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}
