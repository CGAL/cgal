#include <CGAL/Simple_cartesian.h>
#include <CGAL/Deform_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <fstream>
#include <map>
#include <boost/property_map/property_map.hpp>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator  	  vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor  	  edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator        edge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<edge_descriptor, std::size_t>     Internal_edge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_edge_map>     Edge_index_map;

// a model for SurfaceModelingWeightCalculator, use precomputed weights stored in a map
struct Weights_from_map
{
  Weights_from_map(std::map<edge_descriptor, double>* weight_map) : weight_map(weight_map)
  { }
  double operator()(edge_descriptor e, Polyhedron& /*P*/) {
    return (*weight_map)[e];
  }
  std::map<edge_descriptor, double>* weight_map;
};

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map, CGAL::ORIGINAL_ARAP, Weights_from_map> Deform_mesh;

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr << "Cannot open  data/plane.off" << std::endl;
    return 1;
  }

  std::map<edge_descriptor, double> weight_map;
  // Store all weights
  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = boost::edges(mesh); eb != ee; ++eb)
  {
    weight_map[*eb] = 1.0; // store your precomputed weights
  }

  Internal_vertex_map vertex_index_map;
  Internal_edge_map   edge_index_map;
  Deform_mesh deform_mesh(mesh, Vertex_index_map(vertex_index_map), Edge_index_map(edge_index_map), 5, 1e-4, Weights_from_map(&weight_map));

  // Deform mesh as you wish
}