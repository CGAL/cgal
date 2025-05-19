#include <CGAL/Polygonal_schema.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Curves_on_surface_topology.h>

using namespace CGAL::Surface_mesh_topology;
typedef Polygonal_schema_with_combinatorial_map<>           PS;

///////////////////////////////////////////////////////////////////////////////
void create_mesh_1(PS& ps)
{
  ps.add_facet("b -a g f e");
  ps.add_facet("-b k i -g");
  ps.add_facet("-f -i -h");
  ps.add_facet("a -e h j");
  ps.add_facet("-k -m -l -j");
  ps.add_facet("c -r o m");
  ps.add_facet("-o -p -n");
  ps.add_facet("-d l n -q");
  ps.add_facet("d -c q p r");
  ps.perforate_facet("-i");
  ps.perforate_facet("-o");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh_2(PS& ps)
{
  ps.add_facet("a b -a -b o");
  ps.add_facet("-o -p");
  ps.add_facet("c d -c -d p");
  ps.perforate_facet("-o");
}

///////////////////////////////////////////////////////////////////////////////
void create_mesh_3(PS& ps)
{
  ps.add_facet("a b -a -b s");
  ps.add_facet("c d -c -d -s");
  ps.perforate_facet("s");
  ps.perforate_facet("-s");
}

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
  PS ps[3];
  create_mesh_1(ps[0]);
  create_mesh_2(ps[1]);
  create_mesh_3(ps[2]);

  std::size_t num_ps = sizeof(ps) / sizeof(PS);
  bool res=true;

  for(std::size_t i = 0; i < num_ps; ++i)
  {
    Curves_on_surface_topology<PS> cst(ps[i]);
    Path_on_surface<PS> p1(ps[i]), p2(ps[i]), p3(ps[i]), p4(ps[i]);
    create_path_1(p1);
    create_path_2(p2);
    create_path_3(p3);
    create_path_4(p4);

    if(!cst.is_homotopic_to_simple_cycle(p1))
    {
      std::cout<<"ERROR simplicity_double_torus_with_holes surface"
               << (i+1) << "/test1: "
               <<"Path p1 should be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }

    if(cst.is_homotopic_to_simple_cycle(p2))
    {
      std::cout<<"ERROR simplicity_double_torus_with_holes surface"
               << (i+1) << "/test2: "
               <<"Path p2 should not be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }

    if(!cst.is_homotopic_to_simple_cycle(p3))
    {
      std::cout<<"ERROR simplicity_double_torus_with_holes surface"
               << (i+1) << "/test3: "
               <<"Path p3 should be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }

    if(cst.is_homotopic_to_simple_cycle(p4))
    {
      std::cout<<"ERROR simplicity_double_torus_with_holes surface"
               << (i+1) << "/test4: "
               <<"Path p4 should not be homotopic to a simple cycle"
               <<std::endl;
      res=false;
    }
  }

  if (res)
  {
    std::cout<<"SUCCESS simplicity_double_torus_with_holes; all tests ok."<<std::endl;
    return EXIT_SUCCESS;
  }

  return EXIT_FAILURE;
}
