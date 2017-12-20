#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <fstream>
#include <map>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/graph_traits_TriMesh_ArrayKernelT.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;


int main()
{
  {
    typedef CGAL::Surface_mesh<typename Epic::Point_3> SM;
    SM sm;

    std::ifstream in("data/elephant.off");
    in >> sm;
    PMP::isotropic_remeshing(faces(sm),
                             0.02,
                             sm);
    std::ofstream out("sm.off");
    out << sm << std::endl;
  }

  {
  typedef CGAL::Polyhedron_3<Epic> P;
  std::map<boost::graph_traits<P>::face_descriptor, std::size_t> fim;
  P p; 
  std::ifstream in("data/elephant.off");
  in >> p;
  PMP::isotropic_remeshing(faces(p),
                             0.02,
                           p,
                           PMP::parameters::face_index_map(boost::make_assoc_property_map(fim)));
    std::ofstream out("p.off");
    out << p << std::endl;
  }

  {
    typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OM;

    OM om;
    OpenMesh::IO::read_mesh(om, "data/elephant.off");
    om.request_face_status();
    om.request_edge_status();
    om.request_vertex_status();
    PMP::isotropic_remeshing(faces(om),
                             0.02,
                             om);
    
    om.garbage_collection();
    OpenMesh::IO::write_mesh(om, "pm.off");
  }

  {
    typedef OpenMesh::TriMesh_ArrayKernelT</* MyTraits*/> OM;

    OM om;
    OpenMesh::IO::read_mesh(om, "data/elephant.off");
    om.request_face_status();
    om.request_edge_status();
    om.request_vertex_status();
    PMP::isotropic_remeshing(faces(om),
                             0.02,
                             om);
    
    om.garbage_collection();
    OpenMesh::IO::write_mesh(om, "tm.off");
  }

  return 0;
}
