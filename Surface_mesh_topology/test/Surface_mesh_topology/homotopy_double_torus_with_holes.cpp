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
  p.push_back_by_index(38); // Its starting dart
  for (int i=0; i<18; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<SM>& p)
{
  p.push_back_by_index(50); // Its starting dart
  for (int i=0; i<18; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
void create_path_3(Path_on_surface<SM>& p)
{
  p.push_back_by_index(51); // Its starting dart
  for (int i=0; i<18; ++i)
  { p.extend_positive_turn(2); } // Extend the path
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  SM sm;
  std::ifstream in("data/2torus-3borders.off");
  if (!in.is_open())
  {
    std::cout<<"ERROR reading file data/2torus-3borders.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  in>>sm;

  Curves_on_surface_topology<SM> cst(sm);
  Path_on_surface<SM> p1(sm), p2(sm), p3(sm);
  create_path_1(p1);
  create_path_2(p2);
  create_path_3(p3);

  bool res=true;

  if(cst.is_contractible(p1))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test1: "
             <<"Path p1 should not be contractible"
             <<std::endl;
    res=false;
  }

  if (cst.is_contractible(p2))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test2: "
             <<"Path p2 should not be contractible"
             <<std::endl;
    res=false;
  }

  if (cst.is_contractible(p3))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test3: "
             <<"Path p3 should not be contractible"
             <<std::endl;
    res=false;
  }

  if (!cst.are_freely_homotopic(p1, p2))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test4: "
             <<"Path p1 should be freely homotopic with path p2."
             <<std::endl;
    res=false;
  }

  if (cst.are_freely_homotopic(p1, p3))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test5: "
             <<"Path p1 should not be freely homotopic with path p3."
             <<std::endl;
    res=false;
  }

  if (cst.are_freely_homotopic(p2, p3))
  {
    std::cout<<"ERROR homotopy_double_torus_with_holes test6: "
             <<"Path p2 should not be freely homotopic with path p3."
             <<std::endl;
    res=false;
  }

  if (res)
  {
    std::cout<<"SUCCESS homotopy_double_torus_with_holes; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
///////////////////////////////////////////////////////////////////////////////
