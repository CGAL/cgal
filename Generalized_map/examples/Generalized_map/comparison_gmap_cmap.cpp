#include <CGAL/Generalized_map.h>
#include <CGAL/Combinatorial_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<3>   GMap_3;
typedef CGAL::Combinatorial_map<3> CMap_3;
typedef GMap_3::Dart_handle        Dart_const_handle;

template<typename Map>
void test(Map& m)
{
  // Create two tetrahedra.
  typename Map::Dart_handle d1 = m.make_combinatorial_tetrahedron();
  typename Map::Dart_handle d2 = m.make_combinatorial_tetrahedron();
  m.template sew<3>(d1, d2);

  // Display the map characteristics.
  m.display_characteristics(std::cout);
  std::cout<<", valid="<<m.is_valid()<<std::endl;

  unsigned int res = 0;
  // Iterate through all the darts of the first tetrahedron.
  for (typename Map::template Dart_of_cell_range<3>::const_iterator
         it(m.template darts_of_cell<3>(d1).begin()),
         itend(m.template darts_of_cell<3>(d1).end());
       it!=itend; ++it)
    ++res;

  std::cout<<"Number of darts of the first tetrahedron: "<<res<<std::endl;

  res = 0;
  // Iterate through all the darts of the face incident to d1.
  for (typename Map::template Dart_of_cell_range<2>::const_iterator
         it(m.template darts_of_cell<2>(d1).begin()),
         itend(m.template darts_of_cell<2>(d1).end());
       it!=itend; ++it)
    ++res;

  std::cout<<"Number of darts of the face incident to d1: "<<res<<std::endl;
}

int main()
{
  CMap_3 cm;
  GMap_3 gm;

  std::cout<<"Combinatorial map:"<<std::endl;
  test(cm);

  std::cout<<"Generalized map:"<<std::endl;
  test(gm);

  return EXIT_SUCCESS;
}
