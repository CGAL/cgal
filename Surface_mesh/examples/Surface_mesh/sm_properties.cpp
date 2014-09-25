#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<typename K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;

int main()
{
  Mesh m;
  vertex_descriptor v0 = m.add_vertex(K::Point_3(0,2,0));
  vertex_descriptor v1 = m.add_vertex(K::Point_3(2,2,0));
  vertex_descriptor v2 = m.add_vertex(K::Point_3(0,0,0));
  vertex_descriptor v3 = m.add_vertex(K::Point_3(2,0,0));
  vertex_descriptor v4 = m.add_vertex(K::Point_3(1,1,0));
  m.add_face(v3, v1, v4);
  m.add_face(v0, v4, v1);
  m.add_face(v0, v2, v4);
  m.add_face(v2, v3, v4);

  // add a Boolean property to all faces and initialize it to false
  m.add_property_map<Face_descriptor,bool>("f:my_property", false);
  
  // give each vertex a name, the default is empty
  CGAL::Property_map<Vertex_descriptor,std::string> name 
    = m.add_property_map<Vertex_descriptor,std::string>("v:name", "noname");

  // add some names to the vertices
  name[v0] = "hello";
  name[v2] = "world";
  
  // retrieve the point property
  CGAL::Property_map<Vertex_descriptor, K::Point_3> location = m.points();
  for(Mesh::Vertex_iterator it = m.vertices_begin(); it != m.vertices_end(); ++it) { 
    std::cout << name[*it] << " @ " << location[*it] << std::endl;
  }
  
  // delete the string property again
  m.remove_property_map(name);

  return 0;
}

