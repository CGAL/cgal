
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#if defined(CGAL_USE_OPENMESH)
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#endif

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

template <typename Mesh>
void
test()
{
  Mesh m;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0),m);

  typedef typename boost::property_map<Mesh, CGAL::dynamic_vertex_property_t<int> >::type VIM;
  VIM vim = get(CGAL::dynamic_vertex_property_t<int>(), m);
  put(vim, *(vertices(m).first), 7812);
  std::cout << get(vim, *(vertices(m).first)) << std::endl;

  {
    typedef typename boost::property_map<Mesh, CGAL::dynamic_halfedge_property_t<int> >::type VIM;
    VIM vim = get(CGAL::dynamic_halfedge_property_t<int>(), m);
    put(vim, *(halfedges(m).first), 7812);
    
    std::cout << get(vim, *(halfedges(m).first)) << std::endl;
  }
  {
    typedef typename boost::property_map<Mesh, CGAL::dynamic_edge_property_t<int> >::type VIM;
    VIM vim = get(CGAL::dynamic_edge_property_t<int>(), m);
    put(vim, *(edges(m).first), 7812);
    
    std::cout << get(vim, *(edges(m).first)) << std::endl;
  }
  {
    typedef typename boost::property_map<Mesh, CGAL::dynamic_face_property_t<int> >::type VIM;
    VIM vim = get(CGAL::dynamic_face_property_t<int>(), m);
    put(vim, *(faces(m).first), 7812);
    
    std::cout << get(vim, *(faces(m).first)) << std::endl;
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

