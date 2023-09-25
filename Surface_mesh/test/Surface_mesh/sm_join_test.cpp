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
  // vc should be the number of in-active vertex indices
  std::cout << "vertex freelist" << std::endl;
  auto unused_vertices = sm.vertex_freelist();
  for (auto vd: unused_vertices)
    std::cout << vd << std::endl;
  assert(vc == unused_vertices.size());

  std::cout << "face freelist" << std::endl;
  auto unused_faces = sm.face_freelist();
  for (auto fd: unused_faces)
    std::cout << fd << std::endl;
  assert(fc == unused_faces.size());

  std::cout << "edge freelist" << std::endl;
  auto unused_edges = sm.edge_freelist();
  for (auto ed: unused_edges)
    std::cout << ed << std::endl;
  assert(ec == unused_edges.size());
}


int main()
{
  Sm sm1, sm2;
  {
    std::ifstream in(CGAL::data_file_path("meshes/cube.off"));
    in >> sm1;
    CGAL::Euler::remove_center_vertex(*(halfedges(sm1).first),sm1);
  }
  freelist(sm1,1,5,6);

  {
    std::ifstream in(CGAL::data_file_path("meshes/cube.off"));
    in >> sm2;
    CGAL::Euler::remove_center_vertex(*(halfedges(sm2).first),sm2);
  }

  freelist(sm1,1,5,6);

  sm1.join(sm2);
  freelist(sm1,2,10,12);

  return 0;
}

