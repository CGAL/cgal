#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;

int main()
{
  Mesh m;
  Mesh::Vertex_index u;
  for(unsigned int i=0; i < 5; ++i){
    Mesh::Vertex_index v = m.add_vertex(K::Point_3(0,0,i+1));
    if(i==2) u=v;
  }

  m.remove_vertex(u);

  std::cout << "After insertion of 5 vertices and removal of the 3. vertex\n"
            << "# vertices  / # vertices + # removed vertices = "
            << m.number_of_vertices()
            << " / " << m.number_of_vertices() + m.number_of_removed_vertices() << std::endl;

  std::cout << "Iterate over vertices\n";
  {
    for(vertex_descriptor vd : m.vertices()){
      std::cout << m.point(vd) << std::endl;
    }
  }

  // The status of being used or removed is stored in a property map
  Mesh::Property_map<Mesh::Vertex_index,bool> removed
    = m.property_map<Mesh::Vertex_index,bool>("v:removed").first;


  std::cout << "\nIterate over vertices and deleted vertices\n"
            << "# vertices / # vertices + # removed vertices = "
            << m.number_of_vertices()
            << " / " << m.number_of_vertices() + m.number_of_removed_vertices() << std::endl;
    {
    unsigned int i = 0, end = m.number_of_vertices() + m.number_of_removed_vertices();
    for( ; i < end; ++i) {
      vertex_descriptor vh(i);
      assert(m.is_removed(vh) == removed[vh]);
      std::cout << m.point(vh) << ((m.is_removed(vh)) ? "  R\n" : "\n");
    }
  }

  m.collect_garbage();

  std::cout << "\nAfter garbage collection\n"
            << "# vertices / # vertices + # removed vertices = "
            << m.number_of_vertices()
            << " / " << m.number_of_vertices() + m.number_of_removed_vertices() << std::endl;

 {
   unsigned int i = 0, end = m.number_of_vertices() + m.number_of_removed_vertices();
    for( ; i < end; ++i) {
      vertex_descriptor vh(i);
      std::cout << m.point(vh) << ((m.is_removed(vh)) ? "  R\n" : "\n");
    }
  }

  return 0;
}
