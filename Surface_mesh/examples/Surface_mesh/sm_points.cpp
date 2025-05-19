#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

int main()
{
  Mesh m;
  vertex_descriptor v0 = m.add_vertex(K::Point_3(0,0,0));
  vertex_descriptor v1 = m.add_vertex(K::Point_3(1,0,0));
  vertex_descriptor v2 = m.add_vertex(K::Point_3(0,1,0));
  vertex_descriptor v3 = m.add_vertex(K::Point_3(0,0,1));

  face_descriptor fd = m.add_face(v0, v1, v2);
  m.add_face(v1, v0, v3);

  // Access the point for a given vertex
  for(vertex_descriptor vd : vertices_around_face(m.halfedge(fd), m)){
    std ::cout << m.point(vd) << std::endl;
  }

  // Access the range of all points of the mesh
  for( const K::Point_3& p : m.points()){
    std::cout << p << std::endl;
  }

  Mesh ch;
  CGAL::convex_hull_3(m.points().begin(), m.points().end(), ch);

  return 0;
}
