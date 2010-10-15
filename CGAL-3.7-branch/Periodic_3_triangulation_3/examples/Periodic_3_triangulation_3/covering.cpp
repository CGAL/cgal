#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_traits_3<K> GT;

typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT> PDT;

typedef PDT::Point                  Point;
typedef PDT::Covering_sheets        Covering_sheets;

int main()
{
  PDT T;

  // Input point grid (27 points)
  for (double x=0. ; x < .9 ; x += 0.33) {
    for (double y=0. ; y < .9 ; y += 0.33) {
      for (double z=0. ; z < .9 ; z += 0.33) {
        T.insert(Point(x,y,z));	  
  } } }

  Covering_sheets cs = T.number_of_sheets();
  std::cout<<"Current covering: "<<cs[0]<<' '<<cs[1]<<' '<<cs[2]<<std::endl;
  if ( T.is_triangulation_in_1_sheet() ) {                                      // = true
    bool is_extensible = T.is_extensible_triangulation_in_1_sheet_h1()
      || T.is_extensible_triangulation_in_1_sheet_h2();                         // = false
    T.convert_to_1_sheeted_covering();
    cs = T.number_of_sheets();
    std::cout<<"Current covering: "<<cs[0]<<' '<<cs[1]<<' '<<cs[2]<<std::endl;
    if ( is_extensible )                                                        // = false
      std::cout<<"It is safe to change the triangulation here."<<std::endl;
    else
      std::cout<<"It is NOT safe to change the triangulation here!"<<std::endl;

    T.convert_to_27_sheeted_covering();
    cs = T.number_of_sheets();
    std::cout<<"Current covering: "<<cs[0]<<' '<<cs[1]<<' '<<cs[2]<<std::endl;
  }

  std::cout<<"It is (again) safe to modify the triangulation."<<std::endl;

  return 0;
}
