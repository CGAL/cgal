#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>
#include <map>
#include <CGAL/property_map.h>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Polyhedron>::halfedge_iterator    halfedge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<halfedge_descriptor, std::size_t>     Internal_hedge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_hedge_map>     Hedge_index_map;

// A model of SurfaceMeshDeformationWeights using a map of pre-computed weights
struct Weights_from_map
{
  typedef Polyhedron Halfedge_graph;
  Weights_from_map(std::map<halfedge_descriptor, double>* weight_map) : weight_map(weight_map)
  { }
  template<class VertexPointMap>
  double operator()(halfedge_descriptor e, Polyhedron& /*P*/, VertexPointMap /*vpm*/) {
    return (*weight_map)[e];
  }
  std::map<halfedge_descriptor, double>* weight_map;
};

typedef CGAL::Surface_mesh_deformation<Polyhedron, Vertex_index_map, Hedge_index_map, CGAL::ORIGINAL_ARAP, Weights_from_map> Surface_mesh_deformation;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open  data/plane.off" << std::endl;
    return 1;
  }

  std::map<halfedge_descriptor, double> weight_map;
  // Store all the weights
  for(halfedge_descriptor h : halfedges(mesh))
  {
    weight_map[h] = 1.0; // store some precomputed weights
  }

  // Create and initialize the vertex index map
  Internal_vertex_map internal_vertex_index_map;
  Vertex_index_map vertex_index_map(internal_vertex_index_map);
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(vertex_descriptor v : vertices(mesh)) {
    put(vertex_index_map, v, counter++);
  }

  // Create and initialize the halfedge index map
  Internal_hedge_map internal_hedge_index_map;
  Hedge_index_map hedge_index_map(internal_hedge_index_map);
  counter = 0;
  for(halfedge_descriptor h : halfedges(mesh)) {
    put(hedge_index_map, h, counter++);
  }
  Surface_mesh_deformation deform_mesh(mesh,
                                       vertex_index_map,
                                       hedge_index_map,
                                       get(CGAL::vertex_point, mesh),
                                       Weights_from_map(&weight_map));

  // Deform mesh as desired
  // .....
}
