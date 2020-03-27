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

  bool res=true;

  if (cst.is_contractible(p1))
  {
    std::cout<<"ERROR homotopy_small_cylinder test1: "
             <<"Path p1 should not be contractible"
             <<std::endl;
    res=false;
  }

  if (cst.is_contractible(p2))
  {
    std::cout<<"ERROR homotopy_small_cylinder test2: "
             <<"Path p2 should not be contractible"
             <<std::endl;
    res=false;
  }

  if (!cst.is_contractible(p3))
  {
    std::cout<<"ERROR homotopy_small_cylinder test3: "
             <<"Path p3 should be contractible"
             <<std::endl;
    res=false;
  }

  if (cst.are_freely_homotopic(p1, p2))
  {
    std::cout<<"ERROR homotopy_big_cylinder test4: "
             <<"Path p1 should not be freely homotopic with path p2."
             <<std::endl;
    res=false;
  }

  if (res)
  {
    std::cout<<"SUCCESS homotopy_small_cylinder; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
///////////////////////////////////////////////////////////////////////////////
