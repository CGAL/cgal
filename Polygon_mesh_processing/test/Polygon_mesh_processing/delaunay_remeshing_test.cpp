#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <vector>
#include <fstream>

using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh  = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data_remeshing/anchor_dense.off";

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.1;

  double fdist = (argc > 3) ? std::stod(std::string(argv[3])) : 0.01;

  std::cout << "Detect features..." << std::endl;

  using EIFMap = boost::property_map<Mesh, CGAL::edge_is_feature_t>::type;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, 45, eif);

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  Mesh outmesh;
  PMP::delaunay_remeshing(mesh,
    outmesh,
    PMP::parameters::protect_constraints(true)
    .mesh_edge_size(target_edge_length)
    .mesh_facet_size(target_edge_length)
    .mesh_facet_distance(fdist)
    .edge_is_constrained_map(eif));

  std::cout << "Remeshing with edge_is_constrained_map done." << std::endl;

  auto vpmap = get_property_map(boost::vertex_point, mesh);
  std::vector<std::vector<Point> > segments;
  for (auto e : edges(mesh))
  {
    if (get(eif, e))
    {
      std::vector<Point> s = {get(vpmap, source(e, mesh)), get(vpmap, target(e, mesh))};
      segments.push_back(s);
    }
  }

  Mesh outmesh2;
  PMP::delaunay_remeshing(mesh,
    outmesh2,
    PMP::parameters::protect_constraints(true)
    .mesh_edge_size(target_edge_length)
    .mesh_facet_size(target_edge_length)
    .mesh_facet_distance(fdist)
    .polyline_constraints(segments));

  std::cout << "Remeshing with polyline_constraints done." << std::endl;

  assert(outmesh.number_of_vertices() == outmesh2.number_of_vertices());
  assert(outmesh.number_of_faces() == outmesh2.number_of_faces());

  return 0;
}
