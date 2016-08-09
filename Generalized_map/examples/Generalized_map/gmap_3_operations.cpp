#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_handle      Dart_handle;

int main()
{
  GMap_3 gm;

  // Create one hexahedron.
  Dart_handle d1 = gm.make_combinatorial_hexahedron();
  CGAL_assertion( gm.is_valid() );
  CGAL_assertion( gm.is_volume_combinatorial_hexahedron(d1) );

  // Add two edges along two opposite facets.
  gm.insert_cell_1_in_cell_2(d1,gm.alpha<0,1,0>(d1));
  CGAL_assertion( gm.is_valid() );

  Dart_handle d2=gm.alpha<2,1,0,1,2>(d1);
  gm.insert_cell_1_in_cell_2(d2,gm.alpha<0,1,0>(d2));
  CGAL_assertion( gm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges of the cube.
  std::vector<Dart_handle> path;
  path.push_back(gm.alpha<1>(d1));
  path.push_back(gm.alpha<0,2,1>(d1)->alpha(0)->alpha(2)->alpha(1));
  path.push_back(gm.alpha<0>(d2)->alpha(0));
  path.push_back(gm.alpha<2,1>(d2)->alpha(2)->alpha(1));

  // TODO Dart_handle d3=gm.insert_cell_2_in_cell_3(path.begin(),path.end());
  CGAL_assertion( gm.is_valid() );

  // Display the m characteristics.
  gm.display_characteristics(std::cout) << ", valid=" <<
    gm.is_valid() << std::endl;

  // We use the removal operations to get back to the initial cube.
  //  gm.remove_cell<2>(d3);
  CGAL_assertion( gm.is_valid() );

  //  gm.remove_cell<1>(d1->alpha(1));
  CGAL_assertion( gm.is_valid() );

  //  gm.remove_cell<1>(d2->alpha(0));
  CGAL_assertion( gm.is_valid() );

  // Display the m characteristics.
  gm.display_characteristics(std::cout) << ", valid="
					<< gm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}
