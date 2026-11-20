#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/boost/graph/border.h>
#include <CGAL/boost/graph/helpers.h>
#include <cassert>
#include <cmath>
#include <iterator>
#include <vector>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
namespace PMP = CGAL::Polygon_mesh_processing;

void add_bad_patches(Mesh& mesh)
{
  const auto p = mesh.add_vertex(K::Point_3(0, 0, 0));
  const auto q = mesh.add_vertex(K::Point_3(0, 0, 1));
  const auto m = mesh.add_vertex(K::Point_3(0.005, 0.005, 0.5));
  const auto a = mesh.add_vertex(K::Point_3(-100, 0, 0.5));
  const auto b = mesh.add_vertex(K::Point_3(0, -100, 0.5));
  mesh.add_face(a, p, m);
  mesh.add_face(a, m, q);
  mesh.add_face(b, m, p);
  mesh.add_face(b, q, m);
}

void add_planar_patch(Mesh& mesh)
{
  const auto a = mesh.add_vertex(K::Point_3(10, 0, 0));
  const auto b = mesh.add_vertex(K::Point_3(12, 0, 0));
  const auto c = mesh.add_vertex(K::Point_3(12, 2, 0));
  const auto d = mesh.add_vertex(K::Point_3(10, 2, 0));
  const auto center = mesh.add_vertex(K::Point_3(11, 1, 0));
  mesh.add_face(a, b, center);
  mesh.add_face(b, c, center);
  mesh.add_face(c, d, center);
  mesh.add_face(d, a, center);
}

int main()
{
  for(bool include_planar_patch : {false, true})
  {
    Mesh input, output;
    add_bad_patches(input);
    if(include_planar_patch)
      add_planar_patch(input);
    assert(CGAL::is_valid_polygon_mesh(input));
    auto input_patches = input.add_property_map<Mesh::Face_index, std::size_t>("f:patch").first;
    auto input_corners = input.add_property_map<Mesh::Vertex_index, std::size_t>("v:corner").first;
    auto constrained = input.add_property_map<Mesh::Edge_index, bool>("e:constrained").first;
    auto patch_ids = output.add_property_map<Mesh::Face_index, std::size_t>("f:patch").first;
    PMP::remesh_planar_patches(input, output,
      CGAL::parameters::cosine_of_maximum_angle(std::cos(1.4 * CGAL_PI / 180))
        .face_patch_map(input_patches).vertex_corner_map(input_corners)
        .edge_is_constrained_map(constrained),
      CGAL::parameters::face_patch_map(patch_ids));
    Mesh almost_output;
    std::size_t nb_corners = 0;
    for(auto v : vertices(input))
      if(get(input_corners, v) != std::size_t(-1))
        nb_corners = (std::max)(nb_corners, get(input_corners, v) + 1);
    const bool all_remeshed = PMP::remesh_almost_planar_patches(input, almost_output,
      include_planar_patch ? 3 : 2, nb_corners, input_patches, input_corners, constrained);
    assert(!all_remeshed);
    assert(CGAL::is_valid_polygon_mesh(almost_output));
    assert(num_faces(almost_output) == num_faces(output));
    assert(CGAL::is_valid_polygon_mesh(output));
    assert(num_faces(output) == (include_planar_patch ? 6 : 4));
    auto components = output.add_property_map<Mesh::Face_index, std::size_t>("f:component").first;
    assert(PMP::connected_components(output, components) == (include_planar_patch ? 2 : 1));
    std::vector<Mesh::Halfedge_index> boundaries;
    CGAL::extract_boundary_cycles(output, std::back_inserter(boundaries));
    assert(boundaries.size() == (include_planar_patch ? 2 : 1));
    for(auto f : faces(output))
      assert(get(patch_ids, f) < (include_planar_patch ? 3 : 2));
  }
}
