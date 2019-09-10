#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_handle Dart_const_handle;

int main()
{
  GMap_3 gm;

  // Create two tetrahedra.
  Dart_const_handle dh1 = gm.make_combinatorial_tetrahedron();
  Dart_const_handle dh2 = gm.make_combinatorial_tetrahedron();

  // Display the generalized map characteristics.
  gm.display_characteristics(std::cout);
  std::cout<<", valid="<<gm.is_valid()<<std::endl;

  unsigned int res = 0;
  // Iterate through all the darts of the first tetrahedron.
  // Note that GMap_3::Dart_of_orbit_range<0,1,2> in 3D is equivalent to
  // GMap_3::Dart_of_cell_range<3>.
  for (GMap_3::Dart_of_orbit_range<0,1,2>::const_iterator
       it(gm.darts_of_orbit<0,1,2>(dh1).begin()),
       itend(gm.darts_of_orbit<0,1,2>(dh1).end()); it!=itend; ++it)
    ++res;

  std::cout<<"Number of darts of the first tetrahedron: "<<res<<std::endl;

  res = 0;
  // Iterate through all the darts of the face incident to dh2.
  for (GMap_3::Dart_of_orbit_range<0,1>::const_iterator
       it(gm.darts_of_orbit<0,1>(dh2).begin()),
       itend(gm.darts_of_orbit<0,1>(dh2).end()); it!=itend; ++it)
    ++res;

  std::cout<<"Number of darts of the face incident to dh2: "<<res<<std::endl;

  return EXIT_SUCCESS;
}

