#define CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Polygon_mesh_processing/repair_manifoldness.h>

#include <cassert>
#include <fstream>
#include <map>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;

typedef CGAL::Surface_mesh<K::Point_3>                            Surface_mesh;
typedef CGAL::Polyhedron_3<K>                                     Polyhedron;

typedef std::vector<std::vector<std::size_t> >                    Vertices_to_merge_container;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename PolygonMesh>
void read_mesh(const char* fname,
               PolygonMesh& mesh)
{
  std::ifstream input(fname);
  if (!input || !(input >> mesh) || mesh.is_empty())
  {
    std::cerr << fname << " is not a valid off file.\n";
    std::exit(1);
  }
}

template <typename PolygonMesh>
void repair_non_manifold_vertices(const char* fname)
{
  std::cout << "Test: " << fname << std::endl;

  PolygonMesh pmesh;
  read_mesh(fname, pmesh);

  PolygonMesh pmesh_cpy = pmesh;

  // Fix manifoldness by splitting non-manifold vertices
//  std::cout << "-- WITH MERGE --" << std::endl;
//  PMP::treat_non_manifold_vertices(pmesh);
//  std::ofstream("results/merged.off") << pmesh;

  std::cout << "-- WITH SEPARATE --" << std::endl;
  PMP::treat_non_manifold_vertices(pmesh_cpy, PMP::SEPARATE);
  std::ofstream("results/separated.off") << pmesh_cpy;
}

template <typename PolygonMesh>
void repair_non_manifold_vertices()
{
  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_adjacency.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_border.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_simple.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_open_star.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_pinched.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/three_triangles_sharing_a_vertex.off");

  repair_non_manifold_vertices<PolygonMesh>("data_repair/nm_vertices_real.off");
  repair_non_manifold_vertices<PolygonMesh>("data_repair/torso.off");
}

int main(int /*argc*/, char** /*argv*/)
{
  std::cout << "Test Non-Manifold Vertex Repair Functions (SM)" << std::endl;
  repair_non_manifold_vertices<Surface_mesh>();

  std::cout << "Test Non-Manifold Vertex Repair Functions (Polyhedron)" << std::endl;
//  repair_non_manifold_vertices<Polyhedron>();

  return EXIT_SUCCESS;
}
