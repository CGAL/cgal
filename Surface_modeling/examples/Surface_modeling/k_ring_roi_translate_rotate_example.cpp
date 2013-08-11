#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
// Halfedge adaptors for Polyhedron_3
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

#include <CGAL/Deform_mesh.h>

#include <fstream>
#include <map>
#include <queue>
#include <boost/property_map/property_map.hpp>

typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator      vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor      edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator        edge_iterator;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator    out_edge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<edge_descriptor, std::size_t>     Internal_edge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_edge_map>     Edge_index_map;

typedef CGAL::Deform_mesh<Polyhedron, Vertex_index_map, Edge_index_map> Deform_mesh;

// extract vertices which are at most k (inclusive) far from vertex v
std::vector<vertex_descriptor> extract_k_ring(const Polyhedron &P, vertex_descriptor v, int k)
{
  std::map<vertex_descriptor, int>  D;
  std::vector<vertex_descriptor>    Q;
  Q.push_back(v); D[v] = 0;
  std::size_t current_index = 0;

  int dist_v;
  while( current_index < Q.size() && (dist_v = D[ Q[current_index] ]) < k ) {
    v = Q[current_index++];

    out_edge_iterator e, e_end;
    for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
    {
      vertex_descriptor new_v = boost::target(*e, P);
      if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        Q.push_back(new_v);
      }
    }
  }
  return Q;
}

int main()
{
  Polyhedron mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.empty() ) {
    std::cerr<< "Cannot open  data/plane.off";
    return 1;
  }

  // index maps should contain unique indices with 0 offset
  Internal_vertex_map internal_vertex_index_map;
  Vertex_index_map vertex_index_map(internal_vertex_index_map);
  vertex_iterator vb, ve;
  std::size_t counter = 0;
  for(boost::tie(vb, ve) = boost::vertices(mesh); vb != ve; ++vb, ++counter) {
    put(vertex_index_map, *vb, counter);
  }

  Internal_edge_map internal_edge_index_map;
  Edge_index_map edge_index_map(internal_edge_index_map);
  counter = 0;
  edge_iterator eb, ee;
  for(boost::tie(eb, ee) = boost::edges(mesh); eb != ee; ++eb, ++counter) {
    put(edge_index_map, *eb, counter);
  }
//// PREPROCESS SECTION ////
  Deform_mesh deform_mesh(mesh, vertex_index_map, edge_index_map);

  // insert region of interest
  boost::tie(vb,ve) = boost::vertices(mesh);

  std::vector<vertex_descriptor> roi_map = extract_k_ring(mesh, *boost::next(vb, 47), 9);
  std::vector<vertex_descriptor> handles_1_map = extract_k_ring(mesh, *boost::next(vb, 39), 1);
  std::vector<vertex_descriptor> handles_2_map = extract_k_ring(mesh, *boost::next(vb, 97), 1);

  deform_mesh.insert_roi(roi_map.begin(), roi_map.end());
  deform_mesh.insert_handle(handles_1_map.begin(), handles_1_map.end());
  deform_mesh.insert_handle(handles_2_map.begin(), handles_2_map.end());

  deform_mesh.preprocess();
//// DEFORM SECTION ////

  deform_mesh.translate(handles_1_map.begin(), handles_1_map.end(), Eigen::Vector3d(0,0,1));
   // overrides any previous call

  Eigen::Quaternion<double> quad(0.92, 0, 0, -0.38);
  Eigen::Vector3d vect(0, 0, 0);

  deform_mesh.rotate(handles_1_map.begin(), handles_1_map.end(), Deform_mesh::Point(0,0,0), quad, vect);
  deform_mesh.rotate(handles_2_map.begin(), handles_2_map.end(), Deform_mesh::Point(0,0,0), quad, vect);

  deform_mesh.deform();

  std::ofstream output("data/deform_1.off");
  output << mesh; // save deformed mesh
  output.close();

  // Note that translate and rotate are not cumulative,
  // they just use original positions (positions at the time of construction) of the handles while calculating target positions
  deform_mesh.translate(handles_1_map.begin(), handles_1_map.end(), Eigen::Vector3d(0,0.30,0));
  deform_mesh.translate(handles_2_map.begin(), handles_2_map.end(), Eigen::Vector3d(0,0.30,0));

  deform_mesh.set_iterations(10);
  deform_mesh.set_tolerance(0.0);
  deform_mesh.deform(); // will iterate 10 times

  output.open("data/deform_2.off");
  output << mesh;
}

