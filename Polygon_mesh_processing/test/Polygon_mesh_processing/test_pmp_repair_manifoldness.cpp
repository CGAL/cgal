#define CGAL_PMP_REPAIR_MANIFOLDNESS_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_DEBUG
#define CGAL_PMP_REMOVE_SELF_INTERSECTION_OUTPUT
#define CGAL_PMP_SMOOTHING_DEBUG
#define CGAL_PMP_REPAIR_MANIFOLDNESS_SEPARATE_WITH_SMOOTH

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
void read_mesh(const std::string fname,
               PolygonMesh& mesh)
{
  if(!CGAL::IO::read_polygon_mesh(fname, mesh) || is_empty(mesh))
  {
    std::cerr << fname << " is not a valid input file.\n";
    std::exit(1);
  }
}

std::string extract_filename(std::string fname)
{
  fname = fname.substr(fname.find_last_of("/\\") + 1);
  auto pos = fname.find_last_of('.');
  return fname.substr(0, pos);
}

template <typename PolygonMesh>
std::size_t count_nm_vertices(const PolygonMesh& pmesh,
                              const bool verbose = false)
{
  std::vector<typename boost::graph_traits<PolygonMesh>::halfedge_descriptor> nm_points;

  PMP::geometrically_non_manifold_vertices(pmesh, std::inserter(nm_points, nm_points.end()),
                                           CGAL::parameters::verbose(verbose));

  return nm_points.size();
}

template <typename PolygonMesh>
void repair_nm_with_treatment(const std::string fname,
                              const PMP::NM_TREATMENT treatment)
{
  PolygonMesh pmesh;
  read_mesh(fname, pmesh);

  const std::size_t nmn = count_nm_vertices(pmesh, true);
  std::cout << nmn << " nm umbrellas in input" << std::endl;

  PMP::repair_non_manifoldness(pmesh, treatment);

  std::stringstream oss;
  if(treatment == PMP::SEPARATE)
    oss << "results/" << extract_filename(fname) << "_separated.off" << std::ends;
  else if(treatment == PMP::CLIP)
    oss << "results/" << extract_filename(fname) << "_clip.off" << std::ends;
  else
    oss << "results/" << extract_filename(fname) << "_merge.off" << std::ends;

  std::cout << "write " << oss.str() << std::endl;
  CGAL::IO::write_polygon_mesh(oss.str().c_str(), pmesh, CGAL::parameters::stream_precision(17));

  assert(count_nm_vertices(pmesh) == 0);
}

template <typename PolygonMesh>
void repair_nm(const std::string fname)
{
  std::cout << "Test file: " << fname << std::endl;

  std::cout << "================== SEPARATION ==================" << std::endl;
  repair_nm_with_treatment<PolygonMesh>(fname, PMP::SEPARATE);

  std::cout << "================== CLIP ==================" << std::endl;
  repair_nm_with_treatment<PolygonMesh>(fname, PMP::CLIP);

  std::cout << "================== MERGE ==================" << std::endl;
//  repair_nm_with_treatment<PolygonMesh>(fname, PMP::MERGE);
}

template <typename PolygonMesh>
void repair_nm()
{
  repair_nm<PolygonMesh>("data_repair/nm_vertices_simple.off");
  repair_nm<PolygonMesh>("data_repair/nm_vertices_adjacency.off");
  repair_nm<PolygonMesh>("data_repair/nm_vertices_border.off");
  repair_nm<PolygonMesh>("data_repair/nm_vertices_open_star.off");
  repair_nm<PolygonMesh>("data_repair/nm_vertices_pinched.off");
  repair_nm<PolygonMesh>("data_repair/three_triangles_sharing_a_vertex.off");
  repair_nm<PolygonMesh>("data_repair/nm_closed_cubes.off");

  repair_nm<PolygonMesh>("data_repair/nm_edges_simple.off");
  repair_nm<PolygonMesh>("data_repair/nm_edges_triple.off");
  repair_nm<PolygonMesh>("data_repair/nm_edges_grid.off");

  repair_nm<PolygonMesh>("data_repair/nm_vertices_real.off");
  repair_nm<PolygonMesh>("data_repair/torso_no_iv.off");
}

int main(int argc, char** argv)
{
  std::cout.precision(17);

  if(argc > 1)
    repair_nm<Surface_mesh>(argv[1]);

  std::cout << "Test Non-Manifold Vertex Repair Functions (SM)" << std::endl;
  repair_nm<Surface_mesh>();

//  std::cout << "Test Non-Manifold Vertex Repair Functions (Polyhedron)" << std::endl;
//  repair_nm<Polyhedron>();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
