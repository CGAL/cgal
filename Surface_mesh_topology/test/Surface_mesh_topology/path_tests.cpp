#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h>

#include <iostream>

///////////////////////////////////////////////////////////////////////////////
typedef CGAL::Combinatorial_map<2>                           CMap;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3;
using namespace CGAL::Surface_mesh_topology;
///////////////////////////////////////////////////////////////////////////////
bool basic_tests()
{
  bool res=true;
  CMap cmap;
  cmap.make_combinatorial_hexahedron();

  Path_on_surface<CMap> p1(cmap);
  p1.initialize_random_starting_dart();
  p1.extend_straight_positive(3);
  if (!p1.is_valid() || p1.length()!=4 || !p1.is_closed() || !p1.is_simple())
  {
    std::cerr<<"path_tests ERROR: !p1.is_valid() || p1.length()!=4 || "
             <<"!p1.is_closed() || !p1.is_simple()."<<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p2(p1);
  if (p1!=p2 || !p1.are_paths_equals(p2))
  {
    std::cerr<<"path_tests ERROR: p1!=p2 || !p1.are_paths_equals(p2)."<<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p3(cmap);
  p3.push_back(p1.front(), !p1.front_flip());
  p3.extend_straight_negative(3);
  if (p3.length()!=4 || !p3.is_closed())
  {
    std::cerr<<"path_tests ERROR: p3.length()!=4 || !p3.is_closed()."
             <<std::endl;
    res=false;
  }

  p3.reverse();
  if (p1!=p3 || !p1.are_paths_equals(p3))
  {
    std::cerr<<"path_tests ERROR: p1!=p3 || !p1.are_paths_equals(p3)."<<std::endl;
    res=false;
  }

  p3.extend_positive_turn(0); // half-turn
  p3.extend_straight_negative(3);
  if (p3.length()!=8 || !p3.is_closed() || p3.is_simple())
  {
    std::cerr<<"path_tests ERROR: p3.length()!=8 || !p3.is_closed() || p3.is_simple()."
             <<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p4(p3);
  p4.reverse(); // Here p3==p4 because the path is symmetric (it does a round trip)
  if (p3!=p4 || !p3.are_paths_equals(p4))
  {
    std::cerr<<"path_tests ERROR: p3!=p4 || !p3.are_paths_equals(p4)."<<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p5(cmap);
  p5.push_back(p1.front(), !p1.front_flip());
  for (int i=0; i<3; ++i)
  { p5.extend_positive_turn(); }
  p5.reverse();
  if (p1!=p5 || !p1.are_paths_equals(p5))
  {
    std::cerr<<"path_tests ERROR: p1!=p5 || !p1.are_paths_equals(p5)."<<std::endl;
    res=false;
  }

  p1.push_around_face(0);

  Path_on_surface<CMap> p6(cmap);
  p6.push_back(p1.front(), p1.front_flip());
  p6.extend_straight_positive(2);
  p6.extend_positive_turn(1);
  p6.extend_straight_positive(2);
  if (!p6.is_valid() || !p6.is_simple() || p1!=p6)
  {
    std::cerr<<"path_tests ERROR: !p6.is_simple() || !p6.is_simple() || p1!=p6."<<std::endl;
    res=false;
  }

  internal::Light_MQ<CMap> lmq(p6.get_map());
  internal::Path_on_surface_with_rle<internal::Light_MQ<CMap> > p7(lmq, p6);
  if (!p7.is_valid() || p7.size_of_list()!=3)
  {
    std::cerr<<"path_tests ERROR: !p7.is_valid() || size_of_list()!=3."<<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p8(p7);
  if (!p8.is_valid() || p6!=p8 || !p6.are_paths_equals(p8))
  {
    std::cerr<<"path_tests ERROR: !p8.is_valid() || p6!=p8 || !p6.are_paths_equals(p8)."<<std::endl;
    res=false;
  }

  Path_on_surface<CMap> p9(cmap);
  p9.push_back(p1.front(), p1.front_flip());
  if (!p9.can_be_pushed(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back()))) // 1st
  {
    std::cerr<<"path_tests ERROR: 1st !p9.can_be_pushed."<<std::endl;
    res=false;
  }
  else
  { p9.push_back(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back())); }

  if (!p9.can_be_pushed(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back()))) // 2nd
  {
    std::cerr<<"path_tests ERROR: 2nd !p9.can_be_pushed."<<std::endl;
    res=false;
  }
  else
  { p9.push_back(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back())); }

  if (!p9.can_be_pushed(cmap.beta<1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back()))) // 3rd
  {
    std::cerr<<"path_tests ERROR: 3rd !p9.can_be_pushed."<<std::endl;
    res=false;
  }
  else
  { p9.push_back(cmap.beta<1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back())); }

  if (!p9.can_be_pushed(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back()))) // 4th
  {
    std::cerr<<"path_tests ERROR: 4th !p9.can_be_pushed."<<std::endl;
    res=false;
  }
  else
  { p9.push_back(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back())); }

  if (!p9.can_be_pushed(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back()))) // 5th
  {
    std::cerr<<"path_tests ERROR: 5th !p9.can_be_pushed."<<std::endl;
    res=false;
  }
  else
  { p9.push_back(cmap.beta<1,2,1>(p9.back_flip()?cmap.template beta<2>(p9.back()):p9.back())); }

  if (!p9.is_valid() || p6!=p9 || !p6.are_paths_equals(p9))
  {
    std::cerr<<"path_tests ERROR: !p9.is_valid() || p6!=p9 || !p6.are_paths_equals(p9)."<<std::endl;
    res=false;
  }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
int main(/*int argc, char** argv*/)
{
  int res=EXIT_SUCCESS;

  if (!basic_tests()) { res=EXIT_FAILURE; }

  if (res==EXIT_SUCCESS)
  { std::cout<<"All tests in path_tests are OK."<<std::endl; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
