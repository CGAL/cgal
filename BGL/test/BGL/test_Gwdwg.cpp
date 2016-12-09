#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Graph_with_descriptor_with_graph.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> SM;
typedef CGAL::Graph_with_descriptor_with_graph<SM> Mesh;

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

  std::cout << "before properties" << std::endl;

  boost::property_map<Mesh,CGAL::vertex_point_t>::type ppm;
  ppm = get(CGAL::vertex_point,mesh);

  put(ppm, * vertices(mesh).first, Point_3(1,2,3));
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    std::cout << get(ppm,vd) << std::endl;
  }

  SM sm2;
  sm2 = sm;
  Mesh mesh2(sm2);
  try {
    if( target( *(halfedges(mesh).first), mesh2) == *(vertices(mesh).first)){
      std::cerr << "We should not get here" << std::endl;
    }
  } catch(...){
    std::cerr << "we catched it" << std::endl;
  }
    std::cout << "done" << std::endl;
  return 0;
}

