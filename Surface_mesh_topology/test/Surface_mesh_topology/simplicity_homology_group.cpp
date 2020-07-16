#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>

using namespace CGAL::Surface_mesh_topology;
typedef Polygonal_schema_with_combinatorial_map<>           PS;

///////////////////////////////////////////////////////////////////////////////
void create_mesh1(PS& ps)
{
  ps.add_facet("a b -a -b c d -c -d");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh1_path1(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a b c c -a -a -a b");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh1_path2(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a c a a c -a");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh1_path3(Path_on_surface<PS>& p)
{
  p.push_back_by_label("d -c a a b a d b -c -c d b");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh2(PS& ps)
{
  ps.add_facet("a b -a -b c d -c -d e f -e -f");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh2_path1(Path_on_surface<PS>& p)
{
  p.push_back_by_label("b -c e -f c -f c b e c");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh2_path2(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a f -e a -e f");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh2_path3(Path_on_surface<PS>& p)
{
  p.push_back_by_label("f a f -b -e -e b a");
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  PS ps1, ps2;
  create_mesh1(ps1);
  create_mesh2(ps2);

  bool res=true;

  Curves_on_surface_topology<PS> cst1(ps1), cst2(ps2);
  Path_on_surface<PS> p1_1(ps1), p2_1(ps1), p3_1(ps1),
                      p1_2(ps2), p2_2(ps2), p3_2(ps2);

  create_mesh1_path1(p1_1);
  create_mesh1_path2(p2_1);
  create_mesh1_path3(p3_1);
  create_mesh2_path1(p1_2);
  create_mesh2_path2(p2_2);
  create_mesh2_path3(p3_2);

  if(cst1.is_homotopic_to_simple_cycle(p1_1))
  {
    std::cout<<"ERROR simplicity_homology_group test1"
             <<"Path p1_1 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst1.is_homotopic_to_simple_cycle(p2_1))
  {
    std::cout<<"ERROR simplicity_homology_group test2"
             <<"Path p2_1 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst1.is_homotopic_to_simple_cycle(p3_1))
  {
    std::cout<<"ERROR simplicity_homology_group test3"
             <<"Path p3_1 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst2.is_homotopic_to_simple_cycle(p1_2))
  {
    std::cout<<"ERROR simplicity_homology_group test4"
             <<"Path p1_2 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst2.is_homotopic_to_simple_cycle(p2_2))
  {
    std::cout<<"ERROR simplicity_homology_group test5"
             <<"Path p2_2 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if(cst2.is_homotopic_to_simple_cycle(p3_2))
  {
    std::cout<<"ERROR simplicity_homology_group test6"
             <<"Path p3_2 should not be homotopic to a simple cycle"
             <<std::endl;
    res=false;
  }

  if (res)
  {
    std::cout<<"SUCCESS simplicity_homology_group; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
