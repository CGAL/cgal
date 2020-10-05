#include <vector>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef Mesh::Vertex_index vertex_descriptor;
typedef Mesh::Face_index face_descriptor;
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
  vertex_descriptor u = m.add_vertex(K::Point_3(0,1,0));
  vertex_descriptor v = m.add_vertex(K::Point_3(0,0,0));
  vertex_descriptor w = m.add_vertex(K::Point_3(1,0,0));
  vertex_descriptor x = m.add_vertex(K::Point_3(1,1,0));

  /* face_descriptor f = */ m.add_face(u,v,w,x);

  {
    std::cout << "all vertices " << std::endl;

    // The vertex iterator type is a nested type of the Vertex_range
    Mesh::Vertex_range::iterator  vb, ve;

    Mesh::Vertex_range r = m.vertices();
    // The iterators can be accessed through the C++ range API
    vb = r.begin();
    ve = r.end();
    // or the boost Range API
    vb = boost::begin(r);
    ve = boost::end(r);

    // or with boost::tie, as the CGAL range derives from std::pair
    for(boost::tie(vb, ve) = m.vertices(); vb != ve; ++vb){
            std::cout << *vb << std::endl;
    }

    // Instead of the classical for loop one can use
    // the boost macro for a range
    for(vertex_descriptor vd : m.vertices()){
      std::cout << vd << std::endl;
    }

    // or the C++11 for loop. Note that there is a ':' and not a ',' as in BOOST_FOREACH
    for(vertex_descriptor vd : m.vertices()){
      std::cout << vd << std::endl;
    }

  }

  return 0;
}
