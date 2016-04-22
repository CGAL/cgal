#include <cstdlib>

#ifndef CGAL_NO_DEPRECATED_CODE

#define CGAL_NO_DEPRECATION_WARNINGS 1

#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Cell_attribute.h>

typedef CGAL::Combinatorial_map<3> CMap;

bool test()
{
  CMap map;
  for ( int i=0; i<15; ++i )
  {
    CGAL::make_edge(map);
    CGAL::make_combinatorial_polygon(map, 6);
    CGAL::make_combinatorial_tetrahedron(map);
    CGAL::make_combinatorial_hexahedron(map);
  }
  
  for ( int i=0; i<20; ++i )
  {
    CMap::Dart_handle d1=map.darts().begin();
    while ( !map.is_free<3>(d1) ) ++d1;
    CMap::Dart_handle d2=map.darts().begin();
    while ( !map.is_sewable<3>(d1, d2) ) ++d2;
    map.sew<3>(d1,d2);

    if (CGAL::is_face_combinatorial_polygon(map, d1, 6))
    {}
    if (CGAL::is_volume_combinatorial_hexahedron(map, d1))
    {}
    if (CGAL::is_volume_combinatorial_tetrahedron(map, d1))
    {}

    if (CGAL::is_removable<CMap, 1>(map, d1))
      CGAL::remove_cell<CMap, 1>(map, d1);
    else if (CGAL::is_contractible<CMap, 1>(map, d1))
      CGAL::contract_cell<CMap, 1>(map, d1, true);

    CMap::Dart_handle d3=CGAL::insert_cell_0_in_cell_1<CMap>(map, d2);
    CGAL::insert_cell_0_in_cell_2<CMap>(map, d2);

    CGAL::insert_dangling_cell_1_in_cell_2<CMap>(map, d2);
      
    if (CGAL::is_insertable_cell_1_in_cell_2<CMap>(map, d2, d3))
      CGAL::insert_cell_1_in_cell_2<CMap>(map, d2, d3);

    std::vector<CMap::Dart_handle> adarts;
    adarts.push_back(d2);
    adarts.push_back(d3);
    adarts.push_back(map.beta<1>(d3));
    
    if (CGAL::is_insertable_cell_2_in_cell_3(map, adarts.begin(), adarts.end()))
      CGAL::insert_cell_2_in_cell_3<CMap>(map, adarts.begin(), adarts.end());
  }
  
  return true;
}

int main()
{
  std::cout<<"Combinatorial map deprecated test (v1)."<<std::flush;

  if ( !test() )
  {
    std::cout<<" Failed."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<" Success."<<std::endl;
  return EXIT_SUCCESS;
}

#else // CGAL_NO_DEPRECATED_CODE

int main()
{
  return EXIT_SUCCESS;
}

#endif // CGAL_NO_DEPRECATED_CODE
