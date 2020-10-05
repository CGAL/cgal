#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_deformation.h>

#include <fstream>
#include <map>


typedef CGAL::Simple_cartesian<double>                                   Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3>  Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<Mesh>::vertex_iterator     vertex_iterator;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor     edge_descriptor;


typedef Eigen::Vector3d                                                Vector3d;

typedef CGAL::Surface_mesh_deformation<Mesh>     Surface_mesh_deformation;

// Collect the vertices which are at distance less or equal to k
// from the vertex v in the graph of vertices connected by the edges of P
std::vector<vertex_descriptor> extract_k_ring(const Mesh &P, vertex_descriptor v, int k)
{
  std::map<vertex_descriptor, int>  D;
  std::vector<vertex_descriptor>    Q;
  Q.push_back(v); D[v] = 0;
  std::size_t current_index = 0;

  int dist_v;
  while( current_index < Q.size() && (dist_v = D[ Q[current_index] ]) < k ) {
    v = Q[current_index++];

    for(edge_descriptor e : out_edges(v, P))
    {
      halfedge_descriptor he = halfedge(e, P);
      vertex_descriptor new_v = target(he, P);
      if(D.insert(std::make_pair(new_v, dist_v + 1)).second) {
        Q.push_back(new_v);
      }
    }
  }
  return Q;
}

int main()
{
  Mesh mesh;
  std::ifstream input("data/plane.off");

  if ( !input || !(input >> mesh) || mesh.is_empty() ) {
    std::cerr<< "Cannot open data/plane.off";
    return 1;
  }

  // Create the deformation object
  Surface_mesh_deformation deform_mesh(mesh);

  // Select and insert the vertices of the region of interest
  vertex_iterator vb, ve;
  boost::tie(vb,ve) = vertices(mesh);
  std::vector<vertex_descriptor> roi = extract_k_ring(mesh, *std::next(vb, 47), 9);
  deform_mesh.insert_roi_vertices(roi.begin(), roi.end());

  // Select and insert the control vertices
  std::vector<vertex_descriptor> cvertices_1 = extract_k_ring(mesh, *std::next(vb, 39), 1);
  std::vector<vertex_descriptor> cvertices_2 = extract_k_ring(mesh, *std::next(vb, 97), 1);
  deform_mesh.insert_control_vertices(cvertices_1.begin(), cvertices_1.end());
  deform_mesh.insert_control_vertices(cvertices_2.begin(), cvertices_2.end());

  // Apply a rotation to the control vertices
  Eigen::Quaternion<double> quad(0.92, 0, 0, -0.38);
  deform_mesh.rotate(cvertices_1.begin(), cvertices_1.end(), Vector3d(0,0,0), quad);
  deform_mesh.rotate(cvertices_2.begin(), cvertices_2.end(), Vector3d(0,0,0), quad);

  deform_mesh.deform();

  // Save the deformed mesh
  std::ofstream output("deform_1.off");
  output << mesh;
  output.close();

  // Restore the positions of the vertices
  deform_mesh.reset();

  // Apply a translation on the original positions of the vertices (reset() was called before)
  deform_mesh.translate(cvertices_1.begin(), cvertices_1.end(), Vector3d(0,0.3,0));
  deform_mesh.translate(cvertices_2.begin(), cvertices_2.end(), Vector3d(0,0.3,0));

  // Call the function deform() with one-time parameters:
  // iterate 10 times and do not use energy based termination criterion
  deform_mesh.set_iterations(10);
  deform_mesh.set_tolerance(0.0);
  deform_mesh.deform();

  // Save the deformed mesh
  output.open("deform_2.off");
  output << mesh;
}

