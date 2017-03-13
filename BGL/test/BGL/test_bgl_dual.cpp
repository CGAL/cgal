#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Dual.h>

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3                Point;
typedef CGAL::Surface_mesh<Point>      Mesh;
typedef CGAL::Dual<Mesh>               Dual;
typedef boost::graph_traits<Dual>::face_descriptor face_descriptor;
typedef boost::graph_traits<Dual>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Dual>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Dual>::out_edge_iterator out_edge_iterator;
typedef boost::graph_traits<Dual>::in_edge_iterator in_edge_iterator;

int main()
{
  std::ifstream in("data/primal.off");
  Mesh primal;
  in >> primal;

  Dual dual(primal);
  face_descriptor fd = *faces(dual).first;
  halfedge_descriptor hd = halfedge(fd,dual);
  assert(face(hd,dual) == fd); 
  halfedge_descriptor nhd = next(hd,dual);
  assert(hd != nhd);
  assert(hd == prev(nhd,dual));
  assert(face(nhd,dual) == fd);
  BOOST_FOREACH(halfedge_descriptor lhd, halfedges_around_face(hd,dual)){
    assert(face(lhd,dual) == fd); 
  }
  
  vertex_descriptor vd = *vertices(dual).first;
  
  assert(target(halfedge(vd,dual),dual) == vd);
  BOOST_FOREACH(halfedge_descriptor lhd, halfedges_around_target(halfedge(vd,dual),dual)){
    assert(target(lhd,dual) == vd); 
  }

  {
    out_edge_iterator b,e;
    boost::tie(b,e) = out_edges(vd,dual);
    std::cerr << vd << " " << source(*b,dual) << std::endl;
  }

  {
    in_edge_iterator b,e;
    boost::tie(b,e) = in_edges(vd,dual);
    std::cerr << vd << " " << source(*b,dual) << std::endl;
  }
  std::cerr << "done"<< std::endl;
  return 0;
}
