#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K; // default kernel for OpenMesh point type
typedef CGAL::Exact_predicates_exact_constructions_kernel     EK; // alternatice kernel we want to use
typedef OpenMesh::PolyMesh_ArrayKernelT< >                    Mesh;

typedef boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<EK::Point_3> >::type Exact_point_map;

struct Exact_vertex_point_map
{
  // typedef for the property map
  typedef boost::property_traits<Exact_point_map>::value_type value_type;
  typedef boost::property_traits<Exact_point_map>::reference reference;
  typedef boost::property_traits<Exact_point_map>::key_type key_type;
  typedef boost::read_write_property_map_tag category;

  // exterior references
  Exact_point_map exact_point_map;
  Mesh* tm_ptr;

  // Converters
  CGAL::Cartesian_converter<K, EK> to_exact;
  CGAL::Cartesian_converter<EK, K> to_input;

  Exact_vertex_point_map()
    : tm_ptr(nullptr)
  {}

  Exact_vertex_point_map(const Exact_point_map& ep, Mesh& tm)
    : exact_point_map(ep)
    , tm_ptr(&tm)
  {
    for (key_type v : vertices(tm))
      put(exact_point_map, v, to_exact(get(boost::vertex_point, tm, v)));
  }

  friend
  reference get(const Exact_vertex_point_map& map, key_type k)
  {
    CGAL_precondition(map.tm_ptr!=nullptr);
    return get(map.exact_point_map, k);
  }

  friend
  void put(const Exact_vertex_point_map& map, key_type k, const EK::Point_3& p)
  {
    CGAL_precondition(map.tm_ptr!=nullptr);
    put(map.exact_point_map, k, p);
    // create the input point from the exact one
    put(boost::vertex_point, *map.tm_ptr, k, map.to_input(p));
  }
};

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh1, mesh2;

  OpenMesh::IO::read_mesh(mesh1, filename1);
  OpenMesh::IO::read_mesh(mesh2, filename2);

  Mesh out;

  // Create exact map properties
  Exact_point_map pm_1  = get(CGAL::dynamic_vertex_property_t<EK::Point_3>(), mesh1);
  Exact_point_map pm_2  = get(CGAL::dynamic_vertex_property_t<EK::Point_3>(), mesh2);
  Exact_point_map pm_out  = get(CGAL::dynamic_vertex_property_t<EK::Point_3>(), out);

  // Create exact vertex point map that will provide the point to another kernel and fill the default map
  Exact_vertex_point_map vpm_1(pm_1, mesh1);
  Exact_vertex_point_map vpm_2(pm_2, mesh2);
  Exact_vertex_point_map vpm_out(pm_out, out);

  bool valid_union = PMP::corefine_and_compute_union(mesh1, mesh2, out,
    CGAL::parameters::vertex_point_map(vpm_1),
    CGAL::parameters::vertex_point_map(vpm_2),
    CGAL::parameters::vertex_point_map(vpm_out));

  if (valid_union)
  {
    std::cout << "Union was successfully computed\n";
    OpenMesh::IO::write_mesh(out, "union.off");
    return 0;
  }
  std::cout << "Union could not be computed\n";
  return 1;
}
