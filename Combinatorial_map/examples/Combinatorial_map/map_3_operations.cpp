#include <CGAL/Combinatorial_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle        Dart_handle;

int main()
{
  CMap_3 cm;

  // Create one combinatorial hexahedron.
  Dart_handle dh1 = cm.make_combinatorial_hexahedron();

  // Add two edges along two opposite facets.
  CGAL_assertion( cm.is_insertable_cell_1_in_cell_2
                        (cm.beta(dh1,1),cm.beta(dh1,0)) );

  cm.insert_cell_1_in_cell_2(cm.beta(dh1,1), cm.beta(dh1,0));
  CGAL_assertion( cm.is_valid() );

  Dart_handle dh2=cm.beta(dh1,2,1,1,2);

  CGAL_assertion( cm.is_insertable_cell_1_in_cell_2
                              (dh2,cm.beta(dh2,1,1)) );

  cm.insert_cell_1_in_cell_2(dh2, cm.beta(dh2,1,1));
  CGAL_assertion( cm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges
  // of the hexahedron.
  std::vector<Dart_handle> path;
  path.push_back(cm.beta(dh1,1));
  path.push_back(cm.beta(dh1,0,2,1));
  path.push_back(cm.beta(dh2,0));
  path.push_back(cm.beta(dh2,2,1));

  CGAL_assertion( (cm.is_insertable_cell_2_in_cell_3
                               (path.begin(),path.end())) );

  Dart_handle dh3=cm.insert_cell_2_in_cell_3(path.begin(),path.end());
  CGAL_assertion( cm.is_valid() );

  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid=" <<
    cm.is_valid() << std::endl;

  // We use the removal operations to get back to the initial hexahedron.
  CGAL_assertion( (cm.is_removable<2>(dh3)) );
  cm.remove_cell<2>(dh3);
  CGAL_assertion( cm.is_valid() );

  CGAL_assertion( (cm.is_removable<1>(cm.beta(dh1,1))) );
  cm.remove_cell<1>(cm.beta(dh1,1));
  CGAL_assertion( cm.is_valid() );

  CGAL_assertion( (cm.is_removable<1>(cm.beta(dh2,0))) );
  cm.remove_cell<1>(cm.beta(dh2,0));
  CGAL_assertion( cm.is_valid() );

  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid="
                                        << cm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}

