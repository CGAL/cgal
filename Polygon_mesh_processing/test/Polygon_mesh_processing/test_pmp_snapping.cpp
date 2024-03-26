// #define CGAL_PMP_SNAP_DEBUG_PP

#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/snap.h>
#include <CGAL/Polygon_mesh_processing/internal/Snapping/helper.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/property_map.h>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

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
typedef CGAL::Surface_mesh<EPECK::Point_3>                            Exact_surface_mesh;

template <typename Kernel, typename Mesh>
void test_1()
{
  typedef typename Kernel::FT                                         FT;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor     halfedge_descriptor;

  Mesh fg_source, fg_target;

  // empty meshes
  std::cout << "Empty meshes tests..." << std::endl;

  std::size_t res = PMP::experimental::snap_border_vertices(fg_source, fg_target);
  assert(res == 0);

  std::ifstream source_input("data_snapping/border_snapping_source.off");
  if(!source_input || !(source_input >> fg_source))
  {
    std::cerr << "Error: cannot open source mesh\n";
    return;
  }

  // one empty mesh
  res = PMP::experimental::snap_border_vertices(fg_source, fg_target);
  std::cout << "res: " << res << " (expected 0)" << std::endl;
  assert(res == 0);

  res = PMP::experimental::snap_border_vertices(fg_target, fg_source);
  std::cout << "res: " << res << " (expected 0)" << std::endl;
  assert(res == 0);

  std::ifstream target_input("data_snapping/border_snapping_target.off");
  if(!target_input || !(target_input >> fg_target))
  {
    std::cerr << "Error: cannot open target mesh\n";
    return;
  }

  Mesh fg_source_cpy = fg_source;

  // this epsilon value is too small, nothing happens
  std::cout << "*********************** EPS = 0.000000001 *************** " << std::endl;

  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map_small(0.000000001);
  res = PMP::experimental::snap_border_vertices(fg_source_cpy, tol_map_small, fg_target, tol_map_small);
  std::cout << "res: " << res << " (expected 0)" << std::endl;
  assert(res == 0);

  std::vector<halfedge_descriptor> source_halfedge_range;
  PMP::internal::vertices_as_halfedges(vertices(fg_source_cpy), fg_source_cpy, std::back_inserter(source_halfedge_range));
  std::vector<halfedge_descriptor> target_halfedge_range;
  PMP::internal::vertices_as_halfedges(vertices(fg_target), fg_target, std::back_inserter(target_halfedge_range));

  res = PMP::experimental::snap_vertices(source_halfedge_range, fg_source_cpy, tol_map_small,
                                         target_halfedge_range, fg_target, tol_map_small);

  std::cout << "res: " << res << " (expected 0)" << std::endl;
  assert(res == 0);

  // this epsilon value is too big; everything gets snapped!
  std::cout << "*********************** EPS = 0.1 *************** " << std::endl;

  fg_source_cpy = fg_source;
  std::vector<halfedge_descriptor> border_vertices;
  PMP::border_halfedges(fg_source_cpy, std::back_inserter(border_vertices));

  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map_big(0.1);
  res = PMP::experimental::snap_vertices(border_vertices, fg_source_cpy, tol_map_big,
                                         target_halfedge_range, fg_target, tol_map_big,
                                         params::geom_traits(Kernel()),
                                         params::do_lock_mesh(true));

  std::cout << "res: " << res << " (expected 154)" << std::endl;
  assert(res == 154);

  // this is a good value of 'epsilon', but not all expected vertices are projected
  // because the sampling of the border of the source mesh is not uniform
  std::cout << "*********************** EPS = 0.001 *************** " << std::endl;

  fg_source_cpy = fg_source;
  border_vertices.clear();
  PMP::border_halfedges(fg_source_cpy, std::back_inserter(border_vertices));

  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map_good(0.001);
  res = PMP::experimental::snap_vertices(border_vertices, fg_source_cpy, tol_map_good,
                                         target_halfedge_range, fg_target, tol_map_good,
                                         params::default_values(), params::do_lock_mesh(true));
  std::cout << "res: " << res << " vertices" << std::endl;
  assert(res == 76);

  // this one automatically computes an epsilon bound at each vertex
  std::cout << "*********************** EPS = LOCALLY COMPUTED *************** " << std::endl;

  fg_source_cpy = fg_source;
  border_vertices.clear();
  PMP::border_halfedges(fg_source_cpy, std::back_inserter(border_vertices));

  res = PMP::experimental::snap_vertices(border_vertices, fg_source_cpy,
                                         target_halfedge_range, fg_target,
                                         params::default_values(), params::do_lock_mesh(true));
  std::cout << "res: " << res << " vertices" << std::endl;
  assert(res == 77);

  std::ofstream full_snap_out("snapped_mesh.off");
  full_snap_out << std::setprecision(17) << fg_source_cpy;
}

template <typename Kernel, typename Mesh>
void test_2()
{
  typedef typename Kernel::FT                                         FT;
  typedef typename boost::graph_traits<Mesh>::vertex_descriptor       vertex_descriptor;

  Mesh fg_source, fg_target;

  std::ifstream source_input("data_snapping/border_snapping_source_2.off");
  if(!source_input || !(source_input >> fg_source))
  {
    std::cerr << "Error: cannot open source mesh\n";
    return;
  }

  std::ifstream target_input("data_snapping/border_snapping_target_2.off");
  if(!target_input || !(target_input >> fg_target))
  {
    std::cerr << "Error: cannot open target mesh\n";
    return;
  }

  // This configuration has three vertices close to a target vertex, it just to make sure that
  // if a target vertex is already occupied, the source vertex will go to the next one that is
  // within tolerance and is available
  CGAL::Constant_property_map<vertex_descriptor, FT> tol_map(0.5);
  std::size_t res = PMP::experimental::snap_border_vertices(fg_source, tol_map, fg_target, tol_map);
  std::cout << "res: " << res << " vertices" << std::endl;
  assert(res == 3);

  std::ofstream("out.off") << fg_source;
}

template <typename K, typename Mesh>
void test()
{
  test_1<K, Mesh>();
  test_2<K, Mesh>();
}

int main(int, char**)
{
  std::cout << std::endl << "TEST EPICK SURFACE MESH" << std::endl;
  test<EPICK, Surface_mesh>();

  std::cout << std::endl << "TEST EPICK POLYHEDRON" << std::endl;
  test<EPICK, Polyhedron>();

  std::cout << std::endl << "TEST EPECK SURFACE MESH" << std::endl;
  test<EPECK, Exact_surface_mesh>();

  std::cout << "TEST EPECK POLYHEDRON" << std::endl;
  test<EPECK, Exact_polyhedron>();

  return EXIT_SUCCESS;
}
