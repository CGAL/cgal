#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef CGAL::Polyhedron_3<Kernel>                                   Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_iterator    halfedge_iterator;

// Define the maps
typedef std::map<vertex_descriptor, std::size_t>                  Vertex_id_map;
typedef std::map<halfedge_descriptor, std::size_t>                 Hedge_id_map;
typedef boost::associative_property_map<Vertex_id_map>           Vertex_id_pmap;
typedef boost::associative_property_map<Hedge_id_map>             Hedge_id_pmap;


typedef CGAL::Surface_mesh_deformation<Polyhedron, Vertex_id_pmap, Hedge_id_pmap> Surface_mesh_deformation;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off";
    return 1;
  }

  // Init the indices of the vertices from 0 to num_vertices(mesh)-1
  Vertex_id_map vertex_index_map;
  std::size_t counter = 0;
  for(vertex_descriptor v : vertices(mesh))
    vertex_index_map[v]=counter++;

  // Init the indices of the halfedges from 0 to 2*num_edges(mesh)-1
  Hedge_id_map hedge_index_map;
  counter = 0;
  for(halfedge_descriptor h : halfedges(mesh))
    hedge_index_map[h]=counter++;

  Surface_mesh_deformation deform_mesh( mesh,
                                        Vertex_id_pmap(vertex_index_map),
                                        Hedge_id_pmap(hedge_index_map) );

  // Now deform mesh as desired
  // .....
}
