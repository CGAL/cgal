#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/draw_lcc_with_paths.h>

#include <CGAL/Polyhedron_3.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <CGAL/Surface_mesh_curve_topology.h>
#include <CGAL/Path_generators.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/boost/graph/io.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Kernel::Point_3  Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;

///////////////////////////////////////////////////////////////////////////////
void test(const LCC_3_cmap& mesh)
{
  CGAL::Surface_mesh_curve_topology<LCC_3_cmap> smct(mesh);

  CGAL::Path_on_surface<LCC_3_cmap> p1(mesh); // A first path
  p1.generate_random_closed_path(10);

  CGAL::Path_on_surface<LCC_3_cmap> p2(mesh); // A second path
  p2.generate_random_closed_path(10);

  bool res1=smct.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")<<" contractible."<<std::endl;
  
  bool res2=smct.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res2?"IS":"IS NOT")<<" homotopic with path p2 (green)."<<std::endl;

#ifdef CGAL_USE_BASIC_VIEWER
    std::vector<CGAL::Path_on_surface<LCC_3_cmap> > paths;
    paths.push_back(p1);
    paths.push_back(p2);
    CGAL::draw(mesh, paths); // Enable only if CGAL was compiled with Qt5 */
#endif // CGAL_USE_BASIC_VIEWER
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string file=(argc==1?"data/3torus-smooth.off":argv[1]);
  {
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, file.c_str()))
    {
      std::cout<<"ERROR reading file "<<file<<" for linear cell complex."<<std::endl;
      exit(EXIT_FAILURE);
    }
    test(lcc);
  }
  
  {
    Polyhedron p;
    if (!CGAL::read_off(file, p))
    {
      std::cout<<"ERROR reading file "<<file<<" for polyhedron."<<std::endl;
      exit(EXIT_FAILURE);
    }

    LCC_3_cmap lcc; lcc.import_from_halfedge_graph(p);
    test(lcc);
  }

  {
    SM sm;
    if (!CGAL::read_off(file, sm))
    {
      std::cout<<"ERROR reading file "<<file<<" for surface mesh."<<std::endl;
      exit(EXIT_FAILURE);
    }
    LCC_3_cmap lcc; lcc.import_from_halfedge_graph(sm);
    test(lcc);
  }
  
  return EXIT_SUCCESS;
}
