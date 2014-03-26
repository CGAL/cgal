#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<4> CMap_4;
typedef CMap_4::Dart_handle Dart_handle;

Dart_handle make_triangle(CMap_4& amap)
{
 Dart_handle d1 = amap.create_dart();
 Dart_handle d2 = amap.create_dart();
 Dart_handle d3 = amap.create_dart();
 amap.link_beta<1>(d1,d2);
 amap.link_beta<1>(d2,d3);
 amap.link_beta<1>(d3,d1);
 return d1;
}
 
Dart_handle make_tetrahedral(CMap_4& amap)
{
  Dart_handle d1 = make_triangle(amap);
  Dart_handle d2 = make_triangle(amap);
  Dart_handle d3 = make_triangle(amap);
  Dart_handle d4 = make_triangle(amap);
  amap.link_beta<2>(d1, d2);
  amap.link_beta<2>(d3, amap.beta(d2,0));
  amap.link_beta<2>(amap.beta(d1,1), amap.beta(d3,0));
  amap.link_beta<2>(d4, amap.beta(d2,1));
  amap.link_beta<2>(amap.beta(d4,0), amap.beta(d3,1));
  amap.link_beta<2>(amap.beta(d4,1), amap.beta(d1,0));
  return d1;
}

int main()
{
  CMap_4 cm;
  Dart_handle d1 = make_tetrahedral(cm);
  Dart_handle d2 = make_tetrahedral(cm);
  
  cm.sew<4>(d1,d2);
  
  cm.display_characteristics(std::cout);
  std::cout<<", valid="<<cm.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}

