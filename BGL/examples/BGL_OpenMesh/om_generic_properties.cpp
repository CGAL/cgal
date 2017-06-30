
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/property_map.h>

#include <iostream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point_3;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Mesh;

typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;

int main(int argc, char** argv )
{
  Mesh mesh;
  CGAL::make_triangle(Point_3(0,0,0),Point_3(1,0,0),Point_3(1,1,0),mesh);
  
  typedef boost::property_map<Mesh, boost::vertex_property_t<int> >::type VIM;
  VIM vim = add(boost::vertex_property_t<int>("index"), mesh); 
  put(vim, *(vertices(mesh).first), 7812);
  std::cout << get(vim, *(vertices(mesh).first)) << std::endl;
  remove(vim,mesh);

  {
    typedef boost::property_map<Mesh, boost::halfedge_property_t<int> >::type VIM;
    VIM vim = add(boost::halfedge_property_t<int>("index"), mesh);
    put(vim, *(halfedges(mesh).first), 7812);
    
    std::cout << get(vim, *(halfedges(mesh).first)) << std::endl;
    remove(vim,mesh);
  }

  {

    typedef boost::property_map<Mesh, boost::edge_property_t<int> >::type VIM;
    VIM vim = add(boost::edge_property_t<int>("index"), mesh);
    put(vim, *(edges(mesh).first), 7812);
    std::cout << get(vim, *(edges(mesh).first)) << std::endl;
    remove(vim,mesh);
  }

  {
    typedef boost::property_map<Mesh, boost::face_property_t<int> >::type VIM;
    VIM vim = add(boost::face_property_t<int>("index"), mesh);
    put(vim, *(faces(mesh).first), 7812);
    
    std::cout << get(vim, *(faces(mesh).first)) << std::endl;
    remove(vim,mesh);
  }
  return 0;
}
