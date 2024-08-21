#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/non_rigid_mesh_registration.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/Point_set_3.h>

#include <CGAL/Timer.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using Point = K::Point_3;
using Vector = K::Vector_3;
using FT = K::FT;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {
  const std::string source_fn = CGAL::data_file_path("meshes/bear.off");
  const std::string target_fn = CGAL::data_file_path("meshes/bear_bis.off");
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
      std::cout << v0 << " " << v1 << std::endl;
      correspondences_mesh.push_back(std::make_pair(*(source.vertices_begin() + v0), *(target.vertices_begin() + v1)));
    }
  }

  auto vnm = target.add_property_map<Mesh::Vertex_index, K::Vector_3>("v:normal");
  if (vnm.second)
    PMP::compute_vertex_normals(target, vnm.first);

  Mesh::Property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>> vrm = source.add_property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>>("v:rotation").first;
  Mesh::Property_map<Mesh::Vertex_index, K::Vector_3> vtm = source.add_property_map<Mesh::Vertex_index, K::Vector_3>("v:translation").first;

  FT w1 = 1.0;
  FT w2 = 2.0;
  FT w3 = 500;

  PMP::non_rigid_mesh_to_mesh_registration(source, target, vtm, vrm, CGAL::parameters::point_to_point_energy(w1).point_to_plane_energy(w2).as_rigid_as_possible_energy(w3).correspondences(correspondences_mesh));
  PMP::apply_non_rigid_transformation(source, vtm);
  CGAL::IO::write_polygon_mesh("bear" + std::to_string(w1) + "_" + std::to_string(w2) + "_" + std::to_string(w3) + ".off", source);

  return EXIT_SUCCESS;
}