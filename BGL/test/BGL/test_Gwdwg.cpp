#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Gwdwg.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> SM;
typedef CGAL::Gwdwg<SM> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;

int main()
{
  SM sm;
  Mesh mesh(sm);
  std::ifstream in("data/cube.off");
  in >> sm;

  assert( num_vertices(mesh) == num_vertices(sm) );

  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    halfedge_descriptor hd = halfedge(vd,mesh);
    std::cout << vd << std::endl;
  }

  BOOST_FOREACH(halfedge_descriptor hd, halfedges(mesh)){
    vertex_descriptor vd = target(hd,mesh);
    halfedge_descriptor oh = opposite(hd,mesh);
    std::cout << hd << std::endl;
  }

  BOOST_FOREACH(edge_descriptor ed, edges(mesh)){
    std::cout << ed << std::endl;
  }

  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
    std::cout << fd << std::endl;
  }

  
  CGAL::Euler::split_edge(*(halfedges(mesh).first), mesh);	
  CGAL::Euler::remove_center_vertex(*(halfedges(mesh).first), mesh);	

  std::cerr << "done" << std::endl;
  return 0;
}
