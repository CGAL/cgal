#include <CGAL/Surface_mesh.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         SM;

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void create_path_1(Path_on_surface<SM>& p)
{
  p.push_back_by_index(15); // Its starting dart
  for (int i=0; i<10; ++i)
  { p.extend_positive_turn(2); } // Extend the path
  p.extend_positive_turn(1);
  p.extend_positive_turn(2);
  p.extend_positive_turn(1);
  for (int i=0; i<10; ++i)
  { p.extend_positive_turn(2); } // Extend the path
  p.extend_positive_turn(1);
  p.extend_positive_turn(2);
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<SM>& p)
{
  p.push_back_by_index(136); // Its starting dart
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
void create_path_3(Path_on_surface<SM>& p)
{
  p.push_back_by_index(15); // Its starting dart
  p.extend_positive_turn(2);
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
  p.extend_positive_turn(3);
  p.extend_positive_turn(2);
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  SM sm;
  std::ifstream in("data/rond_point_saucisse.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/rond_point_saucisse.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Curves_on_surface_topology<SM> cst(sm);
  Path_on_surface<SM> p1(sm), p2(sm), p3(sm);
  create_path_1(p1);
  create_path_2(p2);
  create_path_3(p3);


  bool res1=cst.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")
           <<" contractible"<<std::endl;
  std::cout<<(res1?"SUCESS":"FAILURE")<<std::endl;

  bool res2=cst.is_contractible(p2);
  std::cout<<"Path p2 (green) "<<(res2?"IS":"IS NOT")
          <<" contractible"<<std::endl;
  std::cout<<(res2?"FAILURE":"SUCESS")<<std::endl;

  bool res3=cst.is_contractible(p3);
  std::cout<<"Path p3 (orange) "<<(res3?"IS":"IS NOT")
          <<" contractible"<<std::endl;
  std::cout<<(res3?"FAILURE":"SUCESS")<<std::endl;

  bool res4=cst.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res4?"IS":"IS NOT")
           <<" freely homotopic with path p2 (green)."<<std::endl;
  std::cout<<(res4?"FAILURE":"SUCESS")<<std::endl;

  bool res5=cst.are_freely_homotopic(p1, p3);
  std::cout<<"Path p1 (pink) "<<(res5?"IS":"IS NOT")
           <<" freely homotopic with path p3 (orange)."<<std::endl;
  std::cout<<(res5?"FAILURE":"SUCESS")<<std::endl;

  bool res6=cst.are_freely_homotopic(p2, p3);
  std::cout<<"Path p2 (green) "<<(res6?"IS":"IS NOT")
           <<" freely homotopic with path p3 (orange)."<<std::endl;
  std::cout<<(res6?"SUCESS":"FAILURE")<<std::endl;

  return res1 && !res2 && !res3 && !res4 && !res5 && res6;
}
///////////////////////////////////////////////////////////////////////////////
