#include <vector>

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<typename K::Point_3> Mesh;
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



  std::vector<vertex_descriptor> vec;
  using namespace boost::assign;
  vec += u, v, w, x;
  face_descriptor f = m.add_face(vec);

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
    BOOST_FOREACH(vertex_descriptor vd, m.vertices()){
      std::cout << vd << std::endl;
    }

    // or the C+11 for loop. Note that there is a ':' and not a ',' as in BOOST_FOREACH 
    for(vertex_descriptor vd : m.vertices()){
      std::cout << vd << std::endl;
    }
    
  }

  {
    std::cout << "vertices around vertex " << v << std::endl;
    CGAL::Vertex_around_target_circulator<Mesh> vbegin(m.halfedge(v),m), done(vbegin);

    do {
      std::cout << *vbegin++ << std::endl;
    } while(vbegin != done);
  }
   
  { 
    std::cout << "vertices around face " << f << std::endl;
    CGAL::Vertex_around_face_iterator<Mesh> vbegin, vend;
    for(boost::tie(vbegin, vend) = vertices_around_face(m.halfedge(f), m);
        vbegin != vend; 
        ++vbegin){
      std::cout << *vbegin << std::endl;
    }
  }


  return 0;
}
