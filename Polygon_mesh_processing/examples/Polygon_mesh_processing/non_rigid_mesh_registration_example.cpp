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
    std::string line;
    while (std::getline(corr, line)) {
      if (line.size() == 0) continue;
      std::istringstream iss(line);

      std::size_t v0, v1;
      if (iss >> v0 >> v1)
        Mesh source, target;
      correspondences_mesh.push_back(std::make_pair(*(source.vertices_begin() + v0), *(target.vertices_begin() + v1)));
    }

    std::cout << correspondences_mesh.size() << " correspondences provided" << std::endl;
  }

  auto vnm = target.add_property_map<Mesh::Vertex_index, K::Vector_3>("v:normal");
  if (vnm.second)
    PMP::compute_vertex_normals(target, vnm.first);

  typename Mesh::Property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>> vrm = source.add_property_map<Mesh::Vertex_index, CGAL::Aff_transformation_3<K>>("v:rotation").first;
  typename Mesh::Property_map<Mesh::Vertex_index, K::Vector_3> vtm = source.add_property_map<Mesh::Vertex_index, K::Vector_3>("v:transformations").first;

/*
  CGAL::Point_set_3<K::Point_3> points;
  points.reserve(target.num_vertices());

  for (auto v : target.vertices())
    points.insert(target.point(v), get(vnm.first, v));

  std::vector<std::pair<Mesh::Vertex_index, std::size_t>> correspondences_pts;
  correspondences_pts.reserve(correspondences_mesh.size());
  for (auto p : correspondences_mesh)
    correspondences_pts.push_back(std::make_pair(p.first, static_cast<std::size_t>(p.second)));

  PMP::non_rigid_mesh_to_points_registration(source, points, vtm, vrm, correspondences_pts, CGAL::parameters::point_to_plane_energy(2.0).point_to_point_energy(0.1).as_rigid_as_possible_energy(10));
  PMP::apply_non_rigid_transformation(source, vtm, vrm);
  CGAL::IO::write_polygon_mesh("mapped.off", source);*/

  FT w1 = 1.0;
  FT w2 = 2.0;
  FT w3 = 5000000;

  PMP::non_rigid_mesh_to_mesh_registration(source, target, vtm, vrm, correspondences_mesh, CGAL::parameters::point_to_point_energy(w1).point_to_plane_energy(w2).as_rigid_as_possible_energy(w3));
  PMP::apply_non_rigid_transformation(source, vtm, vrm);
  CGAL::IO::write_polygon_mesh("bear" + std::to_string(w1) + "_" + std::to_string(w2) + "_" + std::to_string(w3) + ".off", source);

  return EXIT_SUCCESS;
}