#define CGAL_LINKED_WITH_3MF 1

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/IO/3MF.h>
#include <CGAL/Surface_mesh/IO/3MF.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;



int main()
{
  std::vector<Mesh> meshes(2);
  Mesh& m0 = meshes[0];
  Mesh& m1 = meshes[1];

  {
    vertex_descriptor v0 = m0.add_vertex(K::Point_3(0,0,0));
    vertex_descriptor v1 = m0.add_vertex(K::Point_3(1,0,0));
    vertex_descriptor v2 = m0.add_vertex(K::Point_3(0,1,0));
    vertex_descriptor v3 = m0.add_vertex(K::Point_3(0,0,1));

    face_descriptor fd = m0.add_face(v0, v1, v2);
    m0.add_face(v1, v0, v3);
  }

  {
    vertex_descriptor v0 = m1.add_vertex(K::Point_3(10,0,0));
    vertex_descriptor v1 = m1.add_vertex(K::Point_3(11,0,0));
    vertex_descriptor v2 = m1.add_vertex(K::Point_3(10,1,0));
    vertex_descriptor v3 = m1.add_vertex(K::Point_3(10,0,1));

    face_descriptor fd = m1.add_face(v0, v1, v2);
    m1.add_face(v1, v0, v3);
  }

  std::vector<std::string> names(2);
  names[0] = "m0";
  names[1] = "m1";

  CGAL::IO::write_3MF("m.3mf",
                      meshes,
                      names);

  meshes.clear();
  CGAL::IO::read_3MF("m.3mf", meshes);

  std::cout << meshes[0] << std::endl;

  return 0;
}
