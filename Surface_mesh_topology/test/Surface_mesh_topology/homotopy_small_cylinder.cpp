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
  p.push_back_by_index(2); // Its starting dart
  for (int i=0; i<2; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<SM>& p)
{
  p.push_back_by_index(10);  // Its starting dart
  for (int i=0; i<5; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  SM sm;
  std::ifstream in("data/cylinder-2-borders-12darts.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/cylinder-2-borders-12darts.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Curves_on_surface_topology<SM> cst(sm);
  Path_on_surface<SM> p1(sm), p2(sm), p3(sm);
  create_path_1(p1);
  create_path_2(p2);

  bool res1=cst.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")
           <<" contractible"<<std::endl;
  std::cout<<(res1?"FAILURE":"SUCCESS")<<std::endl;

  bool res2=cst.is_contractible(p2);
  std::cout<<"Path p2 (green) "<<(res2?"IS":"IS NOT")
          <<" contractible"<<std::endl;
  std::cout<<(res2?"FAILURE":"SUCCESS")<<std::endl;

  bool res3=cst.is_contractible(p3);
  std::cout<<"Path p3 (orange) "<<(res3?"IS":"IS NOT")
          <<" contractible"<<std::endl;
  std::cout<<(res3?"SUCESS":"FAILURE")<<std::endl;

  bool res4=cst.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res4?"IS":"IS NOT")
          <<" freely_homotopic with p2 (green)"<<std::endl;
  std::cout<<(res4?"FAILURE":"SUCCESS")<<std::endl;

  return (!res1 && !res2 && res3 && !res4);
}
///////////////////////////////////////////////////////////////////////////////
