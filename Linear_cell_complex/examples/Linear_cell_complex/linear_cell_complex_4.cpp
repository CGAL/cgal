#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <iostream>
#include <vector>
#include <cassert>

typedef CGAL::Linear_cell_complex_for_generalized_map<4,5> LCC_4;
typedef LCC_4::Dart_descriptor Dart_descriptor;
typedef LCC_4::Point           Point;
typedef LCC_4::Vector          Vector;
typedef LCC_4::FT              FT;

int main()
{
  LCC_4 lcc;

  // Create two tetrahedra.
  FT p1[5]={ 0, 0, 0, 0, 0}; std::vector<FT> v1(p1,p1+5);
  FT p2[5]={ 0, 2, 0, 0, 0}; std::vector<FT> v2(p2,p2+5);
  FT p3[5]={ 0, 1, 2, 0, 0}; std::vector<FT> v3(p3,p3+5);
  FT p4[5]={ 2, 1, 0, 0, 0}; std::vector<FT> v4(p4,p4+5);
  FT p5[5]={-1, 0, 0, 0, 0}; std::vector<FT> v5(p5,p5+5);
  FT p6[5]={-1, 2, 0, 0, 0}; std::vector<FT> v6(p6,p6+5);
  FT p7[5]={-1, 1, 2, 0, 0}; std::vector<FT> v7(p7,p7+5);
  FT p8[5]={-3, 1, 2, 0, 0}; std::vector<FT> v8(p8,p8+5);

  Dart_descriptor d1 = lcc.make_tetrahedron(Point(5, v1.begin(), v1.end()),
                                        Point(5, v2.begin(), v2.end()),
                                        Point(5, v3.begin(), v3.end()),
                                        Point(5, v4.begin(), v4.end()));

  Dart_descriptor d2 = lcc.make_tetrahedron(Point(5, v5.begin(), v5.end()),
                                        Point(5, v6.begin(), v6.end()),
                                        Point(5, v7.begin(), v7.end()),
                                        Point(5, v8.begin(), v8.end()));

  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  lcc.sew<4>(d1,d2);

  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  // Add one vertex on the middle of the edge containing dart d1.
  Dart_descriptor d3 = lcc.insert_barycenter_in_cell<1>(d1);
  assert( lcc.is_valid() );

  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  // Add one edge to cut the face containing dart d3 in two.
  Dart_descriptor d4 = lcc.insert_cell_1_in_cell_2(d3, lcc.alpha(d1, 1, 0, 1));
  assert( lcc.is_valid() );

  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  // We use removal operations to get back to the initial configuration.
  lcc.remove_cell<1>(d4);
  assert( lcc.is_valid() );

  lcc.remove_cell<0>(d3);
  assert( lcc.is_valid() );

  lcc.unsew<4>(d1);

  lcc.display_characteristics(std::cout);
  std::cout<<", valid="<<lcc.is_valid()<<std::endl;

  return EXIT_SUCCESS;
}
