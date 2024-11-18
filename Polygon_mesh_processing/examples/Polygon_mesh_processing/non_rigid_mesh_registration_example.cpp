#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/non_rigid_mesh_registration.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Point_set_3.h>

#include <CGAL/Timer.h>
#include <sstream>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using Point = K::Point_3;
using Vector = K::Vector_3;
using FT = K::FT;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int, char**) {
  const std::string source_fn = "data/bear_simple.off";
  const std::string target_fn = "data/bear_bis_simple.off";
  const std::string corr_fn = "data/bear_bear_bis.corr";

  Mesh source, target;
  if (!PMP::IO::read_polygon_mesh(source_fn, source))
  {
    std::cerr << "Invalid input " << source_fn << std::endl;
    return 1;
  }

  if (!PMP::IO::read_polygon_mesh(target_fn, target))
  {
    std::cerr << "Invalid input " << target_fn << std::endl;
    return 1;
  }

  std::vector<std::pair<Mesh::Vertex_index, Mesh::Vertex_index>> correspondences_mesh;
  std::ifstream corr(corr_fn);
  if (corr.is_open()) {
    std::size_t v0, v1;
    while (corr >> v0 >> v1) {
      correspondences_mesh.push_back(std::make_pair(*(source.vertices_begin() + v0), *(target.vertices_begin() + v1)));
    }
  }

  Mesh::Property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>> vrm = source.add_property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>>("v:rotation").first;
  Mesh::Property_map<Mesh::Vertex_index, K::Vector_3> vtm = source.add_property_map<Mesh::Vertex_index, K::Vector_3>("v:translation").first;

  FT w1 = 5;
  FT w2 = 20;
  FT w3 = 500;
  new_arap = true;

  std::ostringstream out;
  out.precision(2);
  if (new_arap)
    out << "bear_" << std::fixed << w1 << "_" << w2 << "_" << w3 << "_new.off";
  else
    out << "bear_" << std::fixed << w1 << "_" << w2 << "_" << w3 << "_old.off";

  std::cout << std::endl << out.str() << std::endl;

  PMP::non_rigid_mesh_to_mesh_registration(source, target, vtm, vrm,
    CGAL::parameters::point_to_point_weight(w1)
    .point_to_plane_weight(w2)
    .as_rigid_as_possible_weight(w3)
    .correspondences(std::cref(correspondences_mesh)));

  PMP::apply_non_rigid_transformation(source, vtm);
  CGAL::IO::write_polygon_mesh(out.str(), source);

  return EXIT_SUCCESS;
}