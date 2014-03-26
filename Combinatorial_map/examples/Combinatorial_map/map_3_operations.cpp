#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle        Dart_handle;

int main()
{
  CMap_3 cm;

  // Create one combinatorial hexahedron.
  Dart_handle dh1 = CGAL::make_combinatorial_hexahedron(cm);

  // Add two edges along two opposite facets.
  CGAL_assertion( CGAL::is_insertable_cell_1_in_cell_2
                        (cm,cm.beta(dh1,1),cm.beta(dh1,0)) );

  CGAL::insert_cell_1_in_cell_2(cm,cm.beta(dh1,1),cm.beta(dh1,0));
  CGAL_assertion( cm.is_valid() );

  Dart_handle dh2=cm.beta(dh1,2,1,1,2);

  CGAL_assertion( CGAL::is_insertable_cell_1_in_cell_2
                        (cm,dh2,cm.beta(dh2,1,1)) );
  
  CGAL::insert_cell_1_in_cell_2(cm,dh2,cm.beta(dh2,1,1));
  CGAL_assertion( cm.is_valid() );

  // Insert a facet along these two new edges plus two initial edges
  // of the hexahedron.
  std::vector<Dart_handle> path;
  path.push_back(cm.beta(dh1,1));
  path.push_back(cm.beta(dh1,0,2,1));
  path.push_back(cm.beta(dh2,0));
  path.push_back(cm.beta(dh2,2,1));
  
  CGAL_assertion( (CGAL::is_insertable_cell_2_in_cell_3
		   (cm,path.begin(),path.end())) );

  Dart_handle dh3=CGAL::insert_cell_2_in_cell_3(cm,path.begin(),path.end());
  CGAL_assertion( cm.is_valid() );
  
  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid=" << 
    cm.is_valid() << std::endl;

  // We use the removal operations to get back to the initial hexahedron.
  CGAL_assertion( (CGAL::is_removable<CMap_3, 2>(cm,dh3)) );
  CGAL::remove_cell<CMap_3,2>(cm,dh3);
  CGAL_assertion( cm.is_valid() );

  CGAL_assertion( (CGAL::is_removable<CMap_3, 1>(cm,cm.beta(dh1,1))) );
  CGAL::remove_cell<CMap_3,1>(cm,cm.beta(dh1,1));
  CGAL_assertion( cm.is_valid() );
  
  CGAL_assertion( (CGAL::is_removable<CMap_3, 1>(cm,cm.beta(dh2,0))) );
  CGAL::remove_cell<CMap_3,1>(cm,cm.beta(dh2,0));
  CGAL_assertion( cm.is_valid() );
  
  // Display the combinatorial map characteristics.
  cm.display_characteristics(std::cout) << ", valid=" 
					<< cm.is_valid() << std::endl;

  return EXIT_SUCCESS;
}

