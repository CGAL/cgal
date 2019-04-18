#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>

#include <CGAL/property_map.h>

#include <iostream>
#include <fstream>
#include <set>
#include <vector>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_exact_constructions_kernel             EPECK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel           EPICK;

typedef CGAL::Polyhedron_3<EPECK, CGAL::Polyhedron_items_with_id_3>   Exact_polyhedron;
typedef CGAL::Polyhedron_3<EPICK, CGAL::Polyhedron_items_with_id_3>   Polyhedron;
typedef CGAL::Surface_mesh<EPICK::Point_3>                            Surface_mesh;

template <typename Kernel, typename Mesh>
void test(const char* filename,
          const double large_tolerance,
          const double small_tolerance)
{
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;
  typedef typename Kernel::FT                                         FT;

  Mesh sm;

  std::ifstream in(filename);
  if(!in || !(in >> sm))
  {
    std::cerr << "Error: cannot open mesh\n";
    return;
  }

  std::cout << "------------------" << std::endl;
  std::cout << "num v/f: " << num_vertices(sm) << " " << num_faces(sm) << std::endl;

  std::size_t res;
  Mesh sm_cpy = sm;

  // zero tolerance, just to test the API
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_zero(0);
  PMP::internal::snap_border_vertices_non_conforming(sm_cpy, sm_cpy);
  PMP::internal::snap_border_vertices_non_conforming(sm_cpy, sm_cpy, tol_pmap_zero);
  PMP::internal::snap_border_vertices_non_conforming(sm_cpy, sm_cpy, tol_pmap_zero,
                                                     CGAL::parameters::geom_traits(Kernel()),
                                                     CGAL::parameters::geom_traits(Kernel()));

  // too big, creates wrong snaps
  sm_cpy = sm;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_large(large_tolerance);
  res = PMP::internal::snap_border_vertices_non_conforming(sm_cpy, tol_pmap_large);
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out1("out1.off");
  out1.precision(17);
  out1 << sm_cpy;
  out1.close();

  // too small
  sm_cpy = sm;
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_pmap_small(small_tolerance);
  res = PMP::internal::snap_border_vertices_non_conforming(sm_cpy, tol_pmap_small,
                                                           CGAL::parameters::geom_traits(Kernel()));
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out2("out2.off");
  out2.precision(17);
  out2 << sm_cpy;
  out2.close();

  // automatically computed, custom tolerance at each vertex
  sm_cpy = sm;
  res = PMP::internal::snap_border_vertices_non_conforming(sm_cpy);
  std::cout << "snapped: " << res << std::endl;

  std::ofstream out3("out3.off");
  out3.precision(17);
  out3 << sm_cpy;
  out3.close();
}

void test(const char* filename,
          const double large_tolerance,
          const double small_tolerance)
{
  std::cout << "######################## TEST FILE: " << filename << " ################## " << std::endl;
  std::cout << "~~~~~~~~~~~ TEST EPECK POLYHEDRON ~~~~~~~~~~~" << std::endl;
  test<EPECK, Exact_polyhedron>(filename, large_tolerance, small_tolerance);

  std::cout << std::endl << "~~~~~~~~~~~ TEST EPICK POLYHEDRON ~~~~~~~~~~~" << std::endl;
  test<EPICK, Polyhedron>(filename, large_tolerance, small_tolerance);

  std::cout << std::endl << "~~~~~~~~~~~ TEST EPICK SURFACE MESH ~~~~~~~~~~~" << std::endl;
  test<EPICK, Surface_mesh>(filename, large_tolerance, small_tolerance);
}

int main(int, char**)
{
  test("data_snapping/non_conform_snapping.off", 0.02, 0.001);
  test("data_snapping/non-conform_snapping-hole.off", 0.02, 0.001);
  test("data_snapping/non_conform_snapping-multiple_ccs.off", 0.02, 0.001);
  test("data_snapping/non-conform_snapping-overlap.off", 0.02, 0.001);

  test("data_snapping/real_data.off", 1., 0.05);
  test("data_snapping/real_data_2.off", 2, 0.1);

  return EXIT_SUCCESS;
}

