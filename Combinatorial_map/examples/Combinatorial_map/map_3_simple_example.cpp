#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_const_handle Dart_const_handle;

  int main() 
  { 
    CMap_3 cm;

    // Create two tetrahedra.  
    Dart_const_handle dh1 = CGAL::make_combinatorial_tetrahedron(cm); 
    Dart_const_handle dh2 = CGAL::make_combinatorial_tetrahedron(cm);
  
    // Display the combinatorial map characteristics.
    cm.display_characteristics(std::cout); 
    std::cout<<", valid="<<cm.is_valid()<<std::endl;
  
    unsigned int res = 0; 
    // Iterate over all the darts of the first tetrahedron.  
    // Note that CMap_3::Dart_of_orbit_range<1,2> in 3D is equivalent to 
    // CMap_3::Dart_of_cell_range<3>.  
    for  (CMap_3::Dart_of_orbit_range<1,2>::const_iterator
          it(cm.darts_of_orbit<1,2>(dh1).begin()),
          itend(cm.darts_of_orbit<1,2>(dh1).end()); it!=itend; ++it) 
       ++res;

    std::cout<<"Number of darts of the first tetrahedron: "
             <<res<<std::endl;

    res = 0; 
    // Iterate over all the darts of the facet containing dh2.  
    for (CMap_3::Dart_of_orbit_range<1>::const_iterator
         it(cm.darts_of_orbit<1>(dh2).begin()),
         itend(cm.darts_of_orbit<1>(dh2).end()); it!=itend; ++it) 
       ++res;
  
    std::cout<<"Number of darts of the facet containing dh2: "
             <<res<<std::endl; 

    return EXIT_SUCCESS;
}

