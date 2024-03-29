#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/surface_Delaunay_remeshing.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <vector>
#include <fstream>

using K     = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = K::Point_3;
using Mesh  = CGAL::Surface_mesh<Point>;
using Polyhedron = CGAL::Polyhedron_3<K>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  std::cout << "Random seed = "
    << CGAL::get_default_random().get_seed() << std::endl;

  std::string filename = (argc > 1) ? std::string(argv[1])
    : CGAL::data_file_path("meshes/anchor_dense.off");

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(filename, mesh) || !CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  double target_edge_length = (argc > 2) ? std::stod(std::string(argv[2])) : 0.1;
  double fdist = (argc > 3) ? std::stod(std::string(argv[3])) : 0.01;
  double sharp_angle = (argc > 4) ? std::stod(std::string(argv[4])) : 45.;

  std::cout << "Detect features..." << std::endl;

  // automatically detect features
  using EIFMap = boost::property_map<Mesh, CGAL::edge_is_feature_t>::type;
  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PMP::detect_sharp_edges(mesh, sharp_angle, eif);

  // collect features in a vector of segments
  auto vpmap = get_property_map(boost::vertex_point, mesh);
  std::vector<std::vector<Point> > segments;
  for (auto e : edges(mesh))
  {
    if (get(eif, e))
    {
      auto ps = get(vpmap, source(halfedge(e, mesh), mesh));
      auto pt = get(vpmap, target(halfedge(e, mesh), mesh));
      std::vector<Point> s = { ps, pt };
      segments.push_back(s);
    }
  }
  std::cout << "Detect features done." << std::endl;

  std::cout << "Start remeshing of " << filename
    << " (" << num_faces(mesh) << " faces)..." << std::endl;

  std::cout << "***********************" << std::endl;
  std::cout << "Remeshing with edge_is_constrained_map..." << std::endl;

  Mesh outmesh1 = PMP::surface_Delaunay_remeshing(mesh,
    CGAL::parameters::protect_constraints(true)
    .mesh_edge_size(target_edge_length)
    .mesh_facet_size(target_edge_length)
    .mesh_facet_distance(fdist)
    .edge_is_constrained_map(eif));

  std::cout << "Remeshing with edge_is_constrained_map done." << std::endl;

  std::cout << "***********************" << std::endl;
  std::cout << "Remeshing with polyline_constraints..." << std::endl;

  Mesh outmesh2 = PMP::surface_Delaunay_remeshing(mesh,
    CGAL::parameters::protect_constraints(true)
    .mesh_edge_size(target_edge_length)
    .mesh_facet_size(target_edge_length)
    .mesh_facet_distance(fdist)
    .polyline_constraints(segments));

  std::cout << "Remeshing with polyline_constraints done." << std::endl;

  // Check
  const std::size_t nbv1 = std::distance(outmesh1.vertices().begin(), outmesh1.vertices().end());
  const std::size_t nbv2 = std::distance(outmesh2.vertices().begin(), outmesh2.vertices().end());
  assert(nbv1 == nbv2);

  const std::size_t nbf1 = std::distance(outmesh1.faces().begin(), outmesh1.faces().end());
  const std::size_t nbf2 = std::distance(outmesh2.faces().begin(), outmesh2.faces().end());
  assert(nbf1 == nbf2);

  std::cout << "***********************" << std::endl;
  std::cout << "Remeshing with features_angle_bound..." << std::endl;
  Polyhedron outmesh3 = PMP::surface_Delaunay_remeshing<Mesh, Polyhedron>(mesh,
    CGAL::parameters::protect_constraints(true)
    .mesh_edge_size(target_edge_length)
    .mesh_facet_size(target_edge_length)
    .mesh_facet_distance(fdist)
    .features_angle_bound(sharp_angle));

  std::cout << "Remeshing with features_angle_bound done." << std::endl;

  const std::size_t nbv3 = std::distance(std::begin(vertices(outmesh3)), std::end(vertices(outmesh3)));
  const std::size_t nbf3 = std::distance(std::begin(faces(outmesh3)), std::end(faces(outmesh3)));
  assert(nbv2 == nbv3);
  assert(nbf2 == nbf3);

  return 0;
}
