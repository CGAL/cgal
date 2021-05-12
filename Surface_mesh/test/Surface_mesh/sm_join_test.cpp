#define CGAL_SURFACE_MESH_TEST_SUITE 1  // so that we can access the freelists

#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/Euler_operations.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Sm;

typedef boost::graph_traits<Sm>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Sm>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<Sm>::edge_descriptor edge_descriptor;
typedef boost::graph_traits<Sm>::face_descriptor face_descriptor;

void
freelist(const Sm& sm, int vc, int fc, int ec)
{
  std::cout << "vertex freelist" << std::endl;
  vertex_descriptor vd = sm.vertex_freelist();
  while(vd != sm.null_vertex()){
    --vc;
    std::cout << vd << std::endl;
    halfedge_descriptor hd = halfedge(vd,sm);
    vd = vertex_descriptor((Sm::size_type)hd);
  }
  assert(vc == 0);

  std::cout << "face freelist" << std::endl;
  face_descriptor fd = sm.face_freelist();
  while(fd != sm.null_face()){
    --fc;
    std::cout << fd << std::endl;
    halfedge_descriptor hd = halfedge(fd,sm);
    fd = face_descriptor((Sm::size_type)hd);
  }
  assert(fc == 0);

  std::cout << "edge freelist" << std::endl;
  edge_descriptor ed = sm.edge_freelist();
  while(ed != sm.null_edge()){
    --ec;
    std::cout << ed << std::endl;
    halfedge_descriptor hd = next(halfedge(ed,sm),sm);
    ed = edge(hd,sm);
  }
  assert(ec == 0);
}


int main()
{
  Sm sm1, sm2;
  {
    std::ifstream in("cube.off");
    in >> sm1;
    CGAL::Euler::remove_center_vertex(*(halfedges(sm1).first),sm1);
  }
  freelist(sm1,1,5,6);

  {
    std::ifstream in("cube.off");
    in >> sm2;
    CGAL::Euler::remove_center_vertex(*(halfedges(sm2).first),sm2);
  }

  freelist(sm1,1,5,6);

  sm1.join(sm2);
  freelist(sm1,2,10,12);

  return 0;
}

