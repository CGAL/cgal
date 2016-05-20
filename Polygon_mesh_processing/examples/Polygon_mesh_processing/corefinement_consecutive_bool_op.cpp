#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef Mesh::Property_map<vertex_descriptor,EK::Point_3> Exact_point_map;
typedef Mesh::Property_map<vertex_descriptor,bool> Exact_point_computed;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

struct Coref_point_map
{
  // typedef for the property map
  typedef boost::property_traits<Exact_point_map>::value_type value_type;
  typedef boost::property_traits<Exact_point_map>::reference reference;
  typedef boost::property_traits<Exact_point_map>::category category;
  typedef boost::property_traits<Exact_point_map>::key_type key_type;

  // exterior references
  Exact_point_computed& exact_point_computed;
  Exact_point_map& exact_point;
  Mesh& mesh;

  // Converters
  CGAL::Cartesian_converter<EK, K> to_exact;
  CGAL::Cartesian_converter<K, EK> to_input;

  Coref_point_map(Exact_point_computed& epc,
                  Exact_point_map& ep,
                  Mesh& m)
    :exact_point_computed(epc), exact_point(ep), m(mesh)
  {}

  friend
  reference get(const Coref_point_map& map, key_type k)
  {
    // create exact point if it does not exist
    if (!map.exact_point_computed[k]){
      map.exact_point[k]=to_exact(k);
      map.exact_point_computed[k]=true;
    }
    return map.exact_point[k];
  }

  friend
  void put(const Coref_point_map& map, key_type k, const EK::Point_3& p)
  {
    map.exact_point_computed[k]=true;
    map.exact_point[k]=p;
    // create the input point from the exact one
    map.mesh.point(k)=to_input(p);
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
    mesh1.add_property_map<vertex_descriptor,EK::Point_3>("e:exact_point");
  Exact_point_computed mesh1_exact_points_computed =
    mesh1.add_property_map<vertex_descriptor,bool>("e:exact_points_computed");

  Exact_point_map mesh2_exact_points =
    mesh2.add_property_map<vertex_descriptor,EK::Point_3>("e:exact_point");
  Exact_point_computed mesh1_exact_points_computed =
    mesh2.add_property_map<vertex_descriptor,bool>("e:exact_points_computed");

  Coref_point_map mesh1_maps(mesh1_exact_points, mesh1_exact_points_computed, mesh1);
  Coref_point_map mesh2_maps(mesh2_exact_points, mesh2_exact_points_computed, mesh2);

  Mesh out;
  if ( PMP::intersection(mesh1,
                         mesh2,
                         mesh1,
                         params::vertex_point_map(mesh1_maps) ) )
  {
    if ( PMP::join(mesh1,
                   mesh2,
                   mesh2,
                   params::vertex_point_map(mesh1_maps),
                   params::vertex_point_map(mesh2_maps) ) )
    {
      std::cout << "Intersection and union were successfully computed\n";
      std::ofstream output("inter_union.off");
      output << mesh2;
      return 0;
    }
    std::cout << "Union could not be computed\n";
    return 1;
  }
  std::cout << "Intersection could not be computed\n";
  return 1;
}
