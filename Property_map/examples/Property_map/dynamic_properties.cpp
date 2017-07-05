
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Mesh;

int main()
{
  Mesh m;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0),m);

  typedef CGAL::dynamic_property_map<Mesh, CGAL::vertex_property_t<int> >::type VIM;
  VIM vim = CGAL::add_property(CGAL::vertex_property_t<int>("index"), m);
  put(vim, *(vertices(m).first), 7812);
  std::cout << get(vim, *(vertices(m).first)) << std::endl;
  CGAL::remove_property(vim, m);

  {
    typedef CGAL::dynamic_property_map<Mesh, CGAL::halfedge_property_t<int> >::type VIM;
    VIM vim = CGAL::add_property(CGAL::halfedge_property_t<int>("index"), m);
    put(vim, *(halfedges(m).first), 7812);
    
    std::cout << get(vim, *(halfedges(m).first)) << std::endl;
    CGAL::remove_property(vim,m);
  }
  {
    typedef CGAL::dynamic_property_map<Mesh, CGAL::edge_property_t<int> >::type VIM;
    VIM vim = CGALadd_property(CGAL::edge_property_t<int>("index"), m);
    put(vim, *(edges(m).first), 7812);
    
    std::cout << get(vim, *(edges(m).first)) << std::endl;
    CGAL::remove_property(vim,m);
  }
  {
    typedef CGAL::dynamic_property_map<Mesh, CGAL::face_property_t<int> >::type VIM;
    VIM vim = CGAL::add_property(CGAL::face_property_t<int>("index"), m);
    put(vim, *(faces(m).first), 7812);
    
    std::cout << get(vim, *(faces(m).first)) << std::endl;
    CGAL::remove_property(vim,m);
  }
  return 0;
}

