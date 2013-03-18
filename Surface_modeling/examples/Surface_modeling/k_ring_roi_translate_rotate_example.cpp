#include <CGAL/Deform_mesh.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Eigen_solver_traits.h>

#include <fstream>
#include <map>
#include <queue>
#include <boost/property_map/property_map.hpp>

#include <Eigen/SuperLUSupport>
  
typedef CGAL::Eigen_solver_traits<Eigen::SuperLU<CGAL::Eigen_sparse_matrix<double>::EigenType> > DefaultSolver;
  
typedef CGAL::Simple_cartesian<double>   Kernel;
typedef CGAL::Polyhedron_3<Kernel>       Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator  	  vertex_iterator; 
typedef boost::graph_traits<Polyhedron>::edge_descriptor  	  edge_descriptor;
typedef boost::graph_traits<Polyhedron>::out_edge_iterator    out_edge_iterator;

typedef std::map<vertex_descriptor, std::size_t>   Internal_vertex_map;
typedef std::map<edge_descriptor, std::size_t>     Internal_edge_map;

typedef boost::associative_property_map<Internal_vertex_map>   Vertex_index_map;
typedef boost::associative_property_map<Internal_edge_map>     Edge_index_map;

typedef CGAL::Deform_mesh<Polyhedron, DefaultSolver, Vertex_index_map, Edge_index_map> Deform_mesh;

// extract vertices which are at most k (inclusive) far from vertex v
std::map<vertex_descriptor, int> extract_k_ring(const Polyhedron &P, vertex_descriptor v, int k)
{
  std::map<vertex_descriptor, int>  D;
  std::queue<vertex_descriptor>     Q;
  Q.push(v); D[v] = 0;

  int dist_v;
  while( !Q.empty() && (dist_v = D[Q.front()]) < k ) {
    v = Q.front();
    Q.pop();

    out_edge_iterator e, e_end;
    for(boost::tie(e, e_end) = boost::out_edges(v, P); e != e_end; e++)
    {
      vertex_descriptor new_v = boost::target(*e, P);
      if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        Q.push(new_v);
      }
    }   
  }
  return D;
}

template<class Iterator>
Iterator next_helper(Iterator it, std::size_t n) { 
  Iterator it_next = it;
  while(n-- > 0) { ++it_next; }
  return it_next;
} 

int main()
{
  Polyhedron mesh;
  std::ifstream("models/plane.off") >> mesh;

  Internal_vertex_map vertex_index_map;
  Internal_edge_map   edge_index_map;
//// PREPROCESS SECTION ////
  Deform_mesh deform_mesh(mesh, Vertex_index_map(vertex_index_map), Edge_index_map(edge_index_map)); 

  // insert region of interest
  vertex_iterator vb, ve;
  boost::tie(vb,ve) = boost::vertices(mesh);

  std::map<vertex_descriptor, int> roi_map = extract_k_ring(mesh, *next_helper(vb, 47), 10);
  std::map<vertex_descriptor, int> handles_1_map = extract_k_ring(mesh, *next_helper(vb, 39), 1);
  std::map<vertex_descriptor, int> handles_2_map = extract_k_ring(mesh, *next_helper(vb, 97), 1);

  for(std::map<vertex_descriptor, int>::iterator it = roi_map.begin(); it != roi_map.end(); ++it) {
    deform_mesh.insert_roi(it->first);
  }

  Deform_mesh::Handle_group handles_1 = deform_mesh.create_handle_group();
  for(std::map<vertex_descriptor, int>::iterator it = handles_1_map.begin(); it != handles_1_map.end(); ++it) {
    deform_mesh.insert_handle(handles_1, it->first);
  }
  
  Deform_mesh::Handle_group handles_2 = deform_mesh.create_handle_group();
  for(std::map<vertex_descriptor, int>::iterator it = handles_2_map.begin(); it != handles_2_map.end(); ++it) {
    deform_mesh.insert_handle(handles_2, it->first);
  }

  deform_mesh.preprocess();
//// DEFORM SECTION ////

  deform_mesh.translate(handles_1, Deform_mesh::Vector(0,0,1)); // any latter calls to translate or rotate will override previous calls
  deform_mesh.translate(handles_1, Deform_mesh::Vector(0,0.30,0)); // overrides any previous call

  deform_mesh.deform();

  std::ofstream("deform_1.off") << mesh; // save deformed mesh

  deform_mesh.translate(handles_2, Deform_mesh::Vector(0,0.30,0));

  deform_mesh.set_iterations(10);
  deform_mesh.set_tolerance(0.0);
  deform_mesh.deform(); // will iterate 10 times

  std::ofstream("deform_2.off") << mesh;
}

