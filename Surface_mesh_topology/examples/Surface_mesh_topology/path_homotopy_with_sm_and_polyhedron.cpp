#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/io.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Face_graph_wrapper.h>
#include <CGAL/draw_face_graph_with_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
typedef Kernel::Point_3                                     Point_3;
typedef CGAL::Surface_mesh<Point_3>                         SM;
using namespace CGAL::Surface_mesh_topology;

static unsigned int seed; // Use the same seed for all the tests

///////////////////////////////////////////////////////////////////////////////
template<class FaceGraph>
void test(const FaceGraph& mesh, bool draw, const char* title)
{
  CGAL::Random random(seed);
  Curves_on_surface_topology<FaceGraph> cst(mesh);

  Path_on_surface<FaceGraph> p1(mesh); // A first path
  p1.generate_random_closed_path(10, random);

  Path_on_surface<FaceGraph> p2(mesh); // A second path
  p2.generate_random_closed_path(10, random);

  bool res1=cst.is_contractible(p1, true);
  std::cout<<"Path p1 "<<(res1?"IS":"IS NOT")<<" contractible."<<std::endl;

  bool res2=cst.are_freely_homotopic(p1, p2, true);
  std::cout<<"Path p1 "<<(res2?"IS":"IS NOT")<<" homotopic with path p2."<<std::endl;

  if (draw)
  { CGAL::draw(mesh, {p1, p2}, title); }
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  std::string file=(argc==1?"data/elephant.off":argv[1]);
  bool draw=(argc>2?std::string(argv[2])=="-draw":false);
  seed=static_cast<unsigned int>(CGAL::get_default_random().get_int(0,INT_MAX));

  {
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, file.c_str()))
    {
      std::cout<<"ERROR reading file "<<file<<" for linear cell complex."<<std::endl;
      exit(EXIT_FAILURE);
    }
    test(lcc, draw, "Linear cell complex for combinatorial map");
  }

  {
    Polyhedron p;
    if (!CGAL::read_off(file, p))
    {
      std::cout<<"ERROR reading file "<<file<<" for polyhedron."<<std::endl;
      exit(EXIT_FAILURE);
    }
    test(p, draw, "Polyhedron");
  }

  {
    SM sm;
    if (!CGAL::read_off(file, sm))
    {
      std::cout<<"ERROR reading file "<<file<<" for surface mesh."<<std::endl;
      exit(EXIT_FAILURE);
    }
    test(sm, draw, "Surface mesh");
  }

  return EXIT_SUCCESS;
}
///////////////////////////////////////////////////////////////////////////////
