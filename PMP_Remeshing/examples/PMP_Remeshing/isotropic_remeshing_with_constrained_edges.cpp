#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

struct Border_pmap
{
  typedef boost::graph_traits<Mesh>::edge_descriptor key_type;
  typedef bool value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  Border_pmap(const Mesh& m)
    : m_mesh(m) {}

  inline friend value_type get(const Border_pmap& bmp, const key_type& e) { return is_border(e, bmp.m_mesh); }
  void put(Border_pmap&, const key_type&, const value_type&) { CGAL_error_msg("This property map is read-only."); }
  const Mesh& m_mesh;
};

int main(int argc, char* argv[])
{
  const std::string filename = (argc > 1)
                             ? argv[1]
                             : CGAL::data_file_path("meshes/corner_tris_with_hole.off");

  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  using EIFMap = boost::property_map<Mesh, CGAL::edge_is_feature_t>::type;
  using PIMap = boost::property_map<Mesh, CGAL::face_patch_id_t<int>>::type;
  using VIMap = boost::property_map<Mesh, CGAL::vertex_incident_patches_t<int>>::type;

  EIFMap eif = get(CGAL::edge_is_feature, mesh);
  PIMap pid = get(CGAL::face_patch_id_t<int>(), mesh);

  std::size_t number_of_patches =
      PMP::sharp_edges_segmentation(mesh, 90, eif, pid);
  std::cout << "The input mesh has been segmented into " << number_of_patches << " patches." << std::endl;

  // Border edges are protected (i.e. kept identical),
  // Feature edges are constrained (i.e. can only be collapsed or split).
  PMP::isotropic_remeshing(faces(mesh),
                          0.1,
                          mesh,
                          PMP::parameters::edge_is_protected_map(Border_pmap(mesh))
                          .face_patch_map(pid)
                          .edge_is_constrained_map(eif)
                          .number_of_iterations(3));
  std::cout << "Remeshing done." << std::endl;

  CGAL::IO::write_polygon_mesh("out.ply", mesh, CGAL::parameters::stream_precision(17));

  return 0;
}
