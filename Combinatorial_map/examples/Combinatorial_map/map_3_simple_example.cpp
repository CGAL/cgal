#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle Dart_const_handle;

int main()
{
  CMap_3 cm;

  // Create two tetrahedra.
  Dart_const_handle d1 = make_combinatorial_tetrahedron(cm);  
  Dart_const_handle d2 = make_combinatorial_tetrahedron(cm);
  
  // Display the map characteristics.
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;
  
  unsigned int res = 0;
  // Iterate through all the darts of the first tetrahedron.
  // Note that CMap_3::Dart_of_beta_range<1,2> is in 3D equivalent to
  // CMap_3::Dart_of_cell_range<3>.
  for (CMap_3::Dart_of_orbit_range<1,2>::const_iterator
	 it(cm.darts_of_orbit<1,2>(d1).begin()),
	 itend(cm.darts_of_orbit<1,2>(d1).end());
       it!=itend; ++it)
    ++res;

  std::cout<<"Number of darts of the first tetrahedron: "<<res<<std::endl;

  res = 0;
  // Iterate through all the darts of the face incident to d1.
  for (CMap_3::Dart_of_orbit_range<1>::const_iterator
	 it(cm.darts_of_orbit<1>(d1).begin()),
	 itend(cm.darts_of_orbit<1>(d1).end());
       it!=itend; ++it)
    ++res;
  
  std::cout<<"Number of darts of the face incident to d1: "<<res<<std::endl;

  return EXIT_SUCCESS;
}

