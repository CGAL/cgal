#include <CGAL/Simple_cartesian.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#if defined(CGAL_USE_OPENMESH)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#endif

#include <CGAL/boost/graph/io.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>                       Kernel;
typedef Kernel::Point_3                                      Point;

typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef CGAL::Surface_mesh<Point>                            SM;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, MyTraits>::type LCC;

#if defined(CGAL_USE_OPENMESH)

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OMesh;

#endif

template<typename Mesh>
void test_bgl_read_write(const char* filename)
{
  Mesh sm;
  std::ifstream in(filename);
  CGAL::read_off(in,sm);
  CGAL::write_off(std::cout, sm);
}

int main(int argc, char** argv)
{
  const char* filename=(argc>1)?argv[1]:"data/prim.off";

  test_bgl_read_write<Polyhedron>(filename);
  test_bgl_read_write<SM>(filename);
  test_bgl_read_write<LCC>(filename);

#ifdef CGAL_USE_OPENMESH
  test_bgl_read_write<OMesh>(filename);
#endif

  return EXIT_SUCCESS;
}
