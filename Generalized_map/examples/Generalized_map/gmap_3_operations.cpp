#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>
#include <cassert>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_descriptor      Dart_descriptor;

int main()
{
  GMap_3 gm;

  // Create one combinatorial hexahedron.
  Dart_descriptor d1 = gm.make_combinatorial_hexahedron();

  // Add two edges along two opposite facets.
  gm.insert_cell_1_in_cell_2(d1,gm.alpha<0,1,0>(d1));
  assert( gm.is_valid() );

  Dart_descriptor d2=gm.alpha<2,1,0,1,2>(d1);
  gm.insert_cell_1_in_cell_2(d2,gm.alpha<0,1,0>(d2));
  assert( gm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges
  // of the hexahedron.
  std::vector<Dart_descriptor> path;
  path.push_back(gm.alpha<1>(d1));
  path.push_back(gm.alpha<1,0,1,2,1>(d1));
  path.push_back(gm.alpha<1,0>(d2));
  path.push_back(gm.alpha<2,1>(d2));

  Dart_descriptor d3=gm.insert_cell_2_in_cell_3(path.begin(),path.end());
  assert( gm.is_valid() );

  // Display the generalized map characteristics.
  gm.display_characteristics(std::cout) << ", valid=" <<
                                           gm.is_valid() << std::endl;

   // We use the removal operations to get back to the initial hexahedron.
  gm.remove_cell<2>(d3);
  assert( gm.is_valid() );

  gm.remove_cell<1>(gm.alpha<1>(d1));
  assert( gm.is_valid() );

  gm.remove_cell<1>(gm.alpha<1>(d2));
  assert( gm.is_valid() );
  assert( gm.is_volume_combinatorial_hexahedron(d1) );

  // Display the generalized map characteristics.
  gm.display_characteristics(std::cout) << ", valid="
                                        << gm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}
