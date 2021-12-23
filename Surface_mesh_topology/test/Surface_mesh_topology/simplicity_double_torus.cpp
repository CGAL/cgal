#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>

using namespace CGAL::Surface_mesh_topology;
typedef Polygonal_schema_with_combinatorial_map<>           PS;

///////////////////////////////////////////////////////////////////////////////
void create_path_1(Path_on_surface<PS>& p)
{
  p.push_back_by_label("b b a b -d -d -d -c -d a");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<PS>& p)
{
  p.push_back_by_label("b b a b d c d d d a");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_3(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a b b b a b b a b b");
  p.push_back_by_label("d c");
  p.push_back_by_label("d c d c d");
  p.push_back_by_label("-c -d");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_4(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a b b b a b b a b b");
  p.push_back_by_label("d c");
  p.push_back_by_label("-d -c -d -c -d");
  p.push_back_by_label("-c -d");
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  PS ps;
  ps.add_facet("a b -a -b c d -c -d");

  Curves_on_surface_topology<PS> cst(ps);
  Path_on_surface<PS> p1(ps), p2(ps), p3(ps), p4(ps);
  create_path_1(p1);
  create_path_2(p2);
  create_path_3(p3);
  create_path_4(p4);

  bool res=true;

  if(!cst.is_homotopic_to_simple_cycle(p1))
  {
    std::cout<<"ERROR simplicity_double_torus test1: "
             <<"Path p1 should be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst.is_homotopic_to_simple_cycle(p2))
  {
    std::cout<<"ERROR simplicity_double_torus test2: "
             <<"Path p2 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(!cst.is_homotopic_to_simple_cycle(p3))
  {
    std::cout<<"ERROR simplicity_double_torus test3: "
             <<"Path p3 should be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst.is_homotopic_to_simple_cycle(p4))
  {
    std::cout<<"ERROR simplicity_double_torus test4: "
             <<"Path p4 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if (res)
  {
    std::cout<<"SUCCESS simplicity_double_torus; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
