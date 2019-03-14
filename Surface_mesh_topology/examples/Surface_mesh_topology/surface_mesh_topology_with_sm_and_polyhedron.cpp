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
template<typename Mesh>
void test(const Mesh& mesh)
{
  CGAL::Surface_mesh_curve_topology<Mesh> smct(mesh);

  CGAL::Path_on_surface<Mesh> p1(mesh); // A first path
  p1.push_back_by_index(14); // Its starting dart
  for (int i=0; i<7; ++i)
  { p1.extend_positive_turn(2); } // Extend the path
  
  CGAL::Path_on_surface<Mesh> p2(mesh); // A second path
  std::vector<std::size_t> indices={202, 206, 335, 317, 322, 69, 62, 414}; // All the indices of the darts
  for (int i=0; i<indices.size(); ++i)
  { p2.push_back_by_index(indices[i]); }

  bool res1=smct.is_contractible(p1);
  std::cout<<"Path p1 (pink) "<<(res1?"IS":"IS NOT")<<" contractible."<<std::endl;
  
  bool res2=smct.are_freely_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res2?"IS":"IS NOT")<<" homotopic with path p2 (green)."<<std::endl;

  bool res3=smct.are_base_point_homotopic(p1, p2);
  std::cout<<"Path p1 (pink) "<<(res3?"IS":"IS NOT")<<" base point homotopic with path p3 (orange)."<<std::endl;
}
///////////////////////////////////////////////////////////////////////////////
int main()
{
  {
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, "data/double-torus-example.off"))
    {
      std::cout<<"ERROR reading file data/double-torus-example.off for linear cell complex."<<std::endl;
      exit(EXIT_FAILURE);
    }
    test(lcc);

#ifdef CGAL_USE_BASIC_VIEWER
    /*  TODO SOMETHING std::vector<CGAL::Path_on_surface<LCC_3_cmap> > paths;
  paths.push_back(p1);
  paths.push_back(p2);  
  paths.push_back(p3);  
  CGAL::draw(lcc, paths); // Enable only if CGAL was compiled with Qt5 */
#endif // CGAL_USE_BASIC_VIEWER
  }
  
  {
    Polyhedron p;
    if (!CGAL::read_off(std::string("data/double-torus-example.off"), p))
    {
      std::cout<<"ERROR reading file data/double-torus-example.off for polyhedron."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }

  {
    SM sm;
    if (!CGAL::read_off(std::string("data/double-torus-example.off"), sm))
    {
      std::cout<<"ERROR reading file data/double-torus-example.off for surface mesh."<<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  
  return EXIT_SUCCESS;
}
