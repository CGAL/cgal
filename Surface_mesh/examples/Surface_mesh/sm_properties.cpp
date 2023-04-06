#include <string>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/generators.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
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


  // give each vertex a name, the default is empty
  Mesh::Property_map<vertex_descriptor,std::string> name;
  bool created;
  boost::tie(name, created) = m.add_property_map<vertex_descriptor,std::string>("v:name","m1");
  assert(created);
  // add some names to the vertices
  name[v0] = "hello";
  name[v2] = "world";

  {
    // You get an existing property, and created will be false
    Mesh::Property_map<vertex_descriptor,std::string> name;
    bool created;
    boost::tie(name, created) = m.add_property_map<vertex_descriptor,std::string>("v:name", "");
    assert(! created);
  }

  //  You can't get a property that does not exist
  Mesh::Property_map<face_descriptor,std::string> gnus;
  bool found;
  boost::tie(gnus, found) = m.property_map<face_descriptor,std::string>("v:gnus");
  assert(! found);

  // retrieve the point property for which exists a convenience function
  Mesh::Property_map<vertex_descriptor, K::Point_3> location = m.points();
  for(vertex_descriptor vd : m.vertices()) {
    std::cout << name[vd] << " @ " << location[vd] << std::endl;
  }

  Mesh m2;
  CGAL::make_triangle(K::Point_3(0,0,1), K::Point_3(1,0,1),K::Point_3(0,1,1), m2);

  m2.add_property_map<vertex_descriptor,std::string>("v:name","m2");
  Mesh::Property_map<vertex_descriptor,int> index;
  index = m2.add_property_map<vertex_descriptor,int>("v:index",-1).first;
  int i = 0;
  for (auto v : vertices(m2)) {
      index[v] = i++;
  }

  std::cout << "properties of m1:" << std::endl;
  std::vector<std::string> props = m.properties<vertex_descriptor>();
  for(std::string p : props){
    std::cout << p << std::endl;
  }

  m.join(m2);
  std::cout << "properties of m1 after join:" << std::endl;
  for(std::string p : m.properties<vertex_descriptor>()){
    std::cout << p << std::endl;
  }

  for (auto v : vertices(m)) {
    std::cout << name[v] << std::endl;
  }
  // delete the string property again
  m.remove_property_map(name);

  return 0;
}
