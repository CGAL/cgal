
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/Euler_operations.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Polyhedron_3<K> Mesh;

int main()
{
  Mesh m;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0),m);

  typedef boost::property_map<Mesh, boost::vertex_property_t<int> >::type VIM;
  VIM vim = add(boost::vertex_property_t<int>("index",245), m);
  std::cout << get(vim, *(vertices(m).first)) << std::endl;
  put(vim, *(vertices(m).first), 7812);
  std::cout << get(vim, *(vertices(m).first)) << std::endl;
  remove(vim,m);

  {
    typedef boost::property_map<Mesh, boost::halfedge_property_t<int> >::type VIM;
    VIM vim = add(boost::halfedge_property_t<int>("index"), m);
    put(vim, *(halfedges(m).first), 7812);
    
    std::cout << get(vim, *(halfedges(m).first)) << std::endl;
    remove(vim,m);
  }
  {
    typedef boost::property_map<Mesh, boost::edge_property_t<int> >::type VIM;
    VIM vim = add(boost::edge_property_t<int>("index"), m);
    put(vim, *(edges(m).first), 7812);
    
    std::cout << get(vim, *(edges(m).first)) << std::endl;
    remove(vim,m);
  }
  {
    typedef boost::property_map<Mesh, boost::face_property_t<int> >::type VIM;
    VIM vim = add(boost::face_property_t<int>("index"), m);
    put(vim, *(faces(m).first), 7812);
    
    std::cout << get(vim, *(faces(m).first)) << std::endl;
    remove(vim,m);
  }
  return 0;
}

