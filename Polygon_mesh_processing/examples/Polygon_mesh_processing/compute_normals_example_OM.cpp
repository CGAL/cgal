#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <iostream>
#include <fstream>
#include <boost/foreach.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 Point;
typedef K::Vector_3 Vector;

typedef OpenMesh::PolyMesh_ArrayKernelT< > Surface_mesh;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

template <typename Mesh, typename K, typename H>
class OM_pmap {
public:
  typedef boost::read_write_property_map_tag category;

  typedef K key_type;
  typedef  typename H::value_type value_type;

  typedef value_type& reference;

  OM_pmap(Mesh& m, H h)
    : mesh(m), h(h)
  {}

  inline friend reference get(const OM_pmap<Mesh,K,H>& pm, key_type k)
  {
    return pm.mesh.property(pm.h,k);
  }

  inline friend void put(const OM_pmap<Mesh,K,H>& pm, key_type k, const value_type& v)
  {
    pm.mesh.property(pm.h,k) = v;
  }

  reference operator[](key_type k)
  {
    return mesh.property(h,k);
  }

  Mesh& mesh;
  H h;
};


int main(int argc, char* argv[])
{
  const char* filename = (argc > 1) ? argv[1] : "data/eight.off";
  std::ifstream input(filename);

  Surface_mesh mesh;
  OpenMesh::IO::read_mesh(mesh, filename);

  OpenMesh::FPropHandleT<Vector> fnormalsProp;
  OpenMesh::VPropHandleT<Vector> vnormalsProp;
  mesh.add_property(fnormalsProp);
  mesh.add_property(vnormalsProp);

  OM_pmap<Surface_mesh, face_descriptor, OpenMesh::FPropHandleT<Vector> > fnormals(mesh,fnormalsProp);
  OM_pmap<Surface_mesh, vertex_descriptor, OpenMesh::VPropHandleT<Vector> >vnormals(mesh,vnormalsProp);

  Vector v(0, 0, 0);
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
      put(vnormals, vd, v); 
    }
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
      put(fnormals, fd, v); 
    }
  CGAL::Polygon_mesh_processing::compute_normals
    (mesh,
     vnormals,
     fnormals,
     CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)).
     geom_traits(K()));

  std::cout << "Face normals :" << std::endl;
  BOOST_FOREACH(face_descriptor fd, faces(mesh)){
    std::cout << fnormals[fd] << std::endl;
  }
  std::cout << "Normals at vertices :" << std::endl;
  BOOST_FOREACH(vertex_descriptor vd, vertices(mesh)){
    std::cout << vnormals[vd] << std::endl;
  }
  return 0;
}
