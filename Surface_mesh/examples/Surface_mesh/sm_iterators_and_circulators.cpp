#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<typename K::Point_3> Mesh;

int main()
{
  Mesh m;

  // u            x
  // +------------+
  // |            |
  // |            |
  // |      f     |
  // |            |
  // |            |
  // +------------+
  // v            w

  // Add the points as vertices
  Mesh::Vertex_descriptor u = m.add_vertex(K::Point_3(0,1,0));
  Mesh::Vertex_descriptor v = m.add_vertex(K::Point_3(0,0,0));
  Mesh::Vertex_descriptor w = m.add_vertex(K::Point_3(1,0,0));
  Mesh::Vertex_descriptor x = m.add_vertex(K::Point_3(1,1,0));



  std::vector<Mesh::Vertex_descriptor> vec;
  using namespace boost::assign;
  vec += u, v, w, x;
  Mesh::Face_descriptor f = m.add_face(vec);

  { 
    std::cout << "all vertices " << std::endl;
    BOOST_FOREACH(Mesh::vertex_descriptor vd, m.vertices()){
      std::cout << vd << std::endl;
    }
  }

  {
    std::cout << "vertices around vertex " << v << std::endl;
    Mesh::Vertex_around_target_circulator vbegin(m.halfedge(v),m), vend(vbegin);

    do {
      std::cout << *vbegin++ << std::endl;
    } while(vbegin != vend);
  }
   
  { 
    std::cout << "vertices around face " << f << std::endl;
    Mesh::Vertex_around_face_iterator vbegin, vend;
    for(boost::tie(vbegin, vend) = m.vertices_around_face(m.halfedge(f));
        vbegin != vend; 
        ++vbegin){
      std::cout << *vbegin << std::endl;
    }
  }


  return 0;
}
