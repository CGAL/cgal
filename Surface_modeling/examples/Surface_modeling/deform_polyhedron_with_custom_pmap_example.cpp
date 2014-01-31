#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Halfedge adaptors for Polyhedron_3
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor        edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator            edge_iterator;

//define the maps
typedef std::map<vertex_descriptor, std::size_t>                  Vertex_id_map;
typedef std::map<edge_descriptor, std::size_t>                      Edge_id_map;
typedef boost::associative_property_map<Vertex_id_map>           Vertex_id_pmap;
typedef boost::associative_property_map<Edge_id_map>               Edge_id_pmap;


typedef CGAL::Deform_mesh<Polyhedron, Vertex_id_pmap, Edge_id_pmap> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off";
    return 1;
  }

  // Init the indices of the vertices from 0 to num_vertices(mesh)
  Vertex_id_map vertex_index_map;
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb, ++counter)
    vertex_index_map[*vb]=counter;

  //Init the indices of the halfedges from 0 to num_edges(mesh)
  Edge_id_map edge_index_map;
  counter = 0;
  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = boost::edges(mesh); eb != ee; ++eb, ++counter)
    edge_index_map[*eb]=counter;

  Deform_mesh deform_mesh( mesh,
                           Vertex_id_pmap(vertex_index_map),
                           Edge_id_pmap(edge_index_map) );

  // Now deform mesh as you wish
  // .....
}
