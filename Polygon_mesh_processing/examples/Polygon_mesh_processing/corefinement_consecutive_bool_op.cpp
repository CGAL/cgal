#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel EK;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Mesh::Property_map<vertex_descriptor,EK::Point_3> Exact_point_map;
typedef Mesh::Property_map<vertex_descriptor,bool> Exact_point_computed;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

struct Exact_vertex_point_map
{
  // typedef for the property map
  typedef boost::property_traits<Exact_point_map>::value_type value_type;
  typedef boost::property_traits<Exact_point_map>::reference reference;
  typedef boost::property_traits<Exact_point_map>::category category;
  typedef boost::property_traits<Exact_point_map>::key_type key_type;

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
    for (Mesh::Vertex_index v : vertices(tm))
      exact_point_map[v]=to_exact(tm.point(v));
  }

  friend
  reference get(const Exact_vertex_point_map& map, key_type k)
  {
    CGAL_precondition(map.tm_ptr!=nullptr);
    return map.exact_point_map[k];
  }

  friend
  void put(const Exact_vertex_point_map& map, key_type k, const EK::Point_3& p)
  {
    CGAL_precondition(map.tm_ptr!=nullptr);
    map.exact_point_map[k]=p;
    // create the input point from the exact one
    map.tm_ptr->point(k)=map.to_input(p);
  }
};


int main(int argc, char* argv[])
{
  const char* filename1 = (argc > 1) ? argv[1] : "data/blobby.off";
  const char* filename2 = (argc > 2) ? argv[2] : "data/eight.off";
  std::ifstream input(filename1);

  Mesh mesh1, mesh2;
  if (!input || !(input >> mesh1))
  {
    std::cerr << "First mesh is not a valid off file." << std::endl;
    return 1;
  }
  input.close();
  input.open(filename2);
  if (!input || !(input >> mesh2))
  {
    std::cerr << "Second mesh is not a valid off file." << std::endl;
    return 1;
  }

  Exact_point_map mesh1_exact_points =
    mesh1.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;

  Exact_point_map mesh2_exact_points =
    mesh2.add_property_map<vertex_descriptor,EK::Point_3>("v:exact_point").first;

  Exact_vertex_point_map mesh1_vpm(mesh1_exact_points, mesh1);
  Exact_vertex_point_map mesh2_vpm(mesh2_exact_points, mesh2);

  if ( PMP::corefine_and_compute_intersection(mesh1,
                                              mesh2,
                                              mesh1,
                                              params::vertex_point_map(mesh1_vpm),
                                              params::vertex_point_map(mesh2_vpm),
                                              params::vertex_point_map(mesh1_vpm) ) )
  {
    if ( PMP::corefine_and_compute_union(mesh1,
                                         mesh2,
                                         mesh2,
                                         params::vertex_point_map(mesh1_vpm),
                                         params::vertex_point_map(mesh2_vpm),
                                         params::vertex_point_map(mesh2_vpm) ) )
    {
      std::cout << "Intersection and union were successfully computed\n";
      std::ofstream output("inter_union.off");
      output.precision(17);
      output << mesh2;
      return 0;
    }
    std::cout << "Union could not be computed\n";
    return 1;
  }
  std::cout << "Intersection could not be computed\n";
  return 1;
}
