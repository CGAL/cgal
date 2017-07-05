
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#if defined(CGAL_USE_OPENMESH)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#endif

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point_3;

template <typename Mesh>
void
test()
{
  Mesh m;
  //CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0),m);

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
    VIM vim = CGAL::add_property(CGAL::edge_property_t<int>("index"), m);
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
}
  
int main()
{
  
  typedef CGAL::Surface_mesh<Point_3> SM;
  test<SM>();

  typedef CGAL::Polyhedron_3<K> P;
  test<P>();

  
#if defined(CGAL_USE_OPENMESH)
  typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OM;
  test<OM>();
#endif
  
  return 0;
}

