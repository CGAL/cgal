#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>

using namespace CGAL::Surface_mesh_topology;
typedef Polygonal_schema_with_combinatorial_map<>           PS;

///////////////////////////////////////////////////////////////////////////////
void create_mesh_1(PS& ps)
{
  ps.add_facet("a b c d");
  ps.add_facet("-b e -d f");
  ps.add_facet("-a -e");
  ps.add_facet("-c -f");
  ps.perforate_facet("-a");
  ps.perforate_facet("-c");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh_2(PS& ps)
{
  ps.add_facet("a b c d");
  ps.add_facet("-b e -d f");
  ps.add_facet("-a -e");
  ps.add_facet("-c -f");
  ps.perforate_facet("-a");
  ps.perforate_facet("-c");
  ps.perforate_facet("a");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh_3(PS& ps)
{
  ps.add_facet("a b c d");
  ps.add_facet("-b e -d f");
  ps.add_facet("-a -e");
  ps.add_facet("-c -f");
  ps.perforate_facet("-a");
  ps.perforate_facet("-c");
  ps.perforate_facet("e");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_1(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a e");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_2(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a b -f d a e");
}

///////////////////////////////////////////////////////////////////////////////
void create_path_3(Path_on_surface<PS>& p)
{
  p.push_back_by_label("a b c d");
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  PS ps[2];
  create_mesh_1(ps[0]);
  create_mesh_2(ps[1]);

  std::size_t num_ps = sizeof(ps) / sizeof(PS);
  bool res=true;

  for(std::size_t i = 0; i < num_ps; ++i)
  {
    Curves_on_surface_topology<PS> cst(ps[i]);
    Path_on_surface<PS> p1(ps[i]), p2(ps[i]), p3(ps[i]);
    create_path_1(p1);
    create_path_2(p2);
    create_path_3(p3);

    if(!cst.is_homotopic_to_simple_cycle(p1))
    {
      std::cout<<"ERROR simplicity_cylinder surface"
               << (i+1) << "/test1: "
               <<"Path p1 should be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }

    if(cst.is_homotopic_to_simple_cycle(p2))
    {
      std::cout<<"ERROR simplicity_cylinder surface"
               << (i+1) << "/test2: "
               <<"Path p2 should not be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }

    if(!cst.is_homotopic_to_simple_cycle(p3))
    {
      std::cout<<"ERROR simplicity_cylinder surface"
               << (i+1) << "/test3: "
               <<"Path p3 should be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }
  }

  if (res)
  {
    std::cout<<"SUCCESS simplicity_cylinder; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
