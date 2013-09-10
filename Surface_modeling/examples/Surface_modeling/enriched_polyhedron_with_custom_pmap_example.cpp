#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Halfedge adaptors for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

#include <fstream>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Polyhedron_items_with_id_3.h>

// Property map using extra field in KeyType
template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;

    reference operator[](key_type key) const { return key->id(); }
};

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron; //enriched polyhedron

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator        edge_iterator;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off";
    return 1;
  }
  // index maps must contain an index unique per vertex starting from 0
  // to the total number of vertices
  Vertex_index_map vertex_index_map;
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb, ++counter) {
    put(vertex_index_map, *vb, counter);
  }

  Edge_index_map edge_index_map;
  counter = 0;
  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = boost::edges(mesh); eb != ee; ++eb, ++counter) {
    put(edge_index_map, *eb, counter);
  }
  Deform_mesh deform_mesh(mesh, vertex_index_map, edge_index_map);

  // Deform mesh as you wish
}

