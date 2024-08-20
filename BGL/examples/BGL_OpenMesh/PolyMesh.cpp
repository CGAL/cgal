#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>
#include <fstream>

namespace CGAL {
namespace IO {

template <typename SM, typename VSelectionPM, typename EFeaturePM>
bool read_OM(std::string fname, SM& sm, VSelectionPM vspm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> Mesh;
  Mesh mesh;
  OpenMesh::IO::Options options  = OpenMesh::IO::Options::Status;
  bool ok = OpenMesh::IO::read_mesh(mesh, fname, options);
  if(! ok){
    return false;
  }

  std::map<vertex_descriptor,sm_vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<halfedge_descriptor,sm_halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  CGAL::copy_face_graph<Mesh,SM>(mesh, sm, CGAL::parameters::vertex_to_vertex_map(v2vpmap).halfedge_to_halfedge_map(h2hpmap));

  if(options.vertex_has_status()){
    for(auto v : vertices(mesh)){
        put(vspm, v2v[v], mesh.status(v).selected());
    }
  }else{
    std::cout << "no vertex status" << std::endl;
  }

  if(options.edge_has_status()){
    for(auto e : edges(mesh)){
        auto sme = edge(h2h[halfedge(e,mesh)], sm);
        put(efpm, sme , mesh.status(e).feature());
    }
  }else{
    std::cout << "no edge status" << std::endl;
  }
  return true;
}

} // namespace IO
} // namespace CGAL


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> Mesh;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<SM>::vertex_descriptor sm_vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SM>::halfedge_descriptor sm_halfedge_descriptor;
typedef boost::graph_traits<SM>::edge_descriptor sm_edge_descriptor;

int main(int argc, char** argv )
{
  Mesh mesh;

  const std::string filename = (argc>1)?argv[1]:CGAL::data_file_path("meshes/in.off");
  const char* outname= (argc>2)?argv[2]:"out.om";
  CGAL::IO::read_polygon_mesh(filename, mesh);

  mesh.request_vertex_status();
  mesh.request_edge_status();

  int i = 0;
  for(auto v : vertices(mesh)){
        mesh.status(v).set_selected((i%2) == 0);
        ++i;
  }

  i = 0;
  for(auto e : edges(mesh)){
        mesh.status(e).set_feature(i > 2);
        ++i;
  }

  OpenMesh::IO::write_mesh(mesh, outname, OpenMesh::IO::Options::Status);

  SM sm;

  std::map<sm_vertex_descriptor,bool> sm_selected_map;
  auto sm_selected_pmap = boost::make_assoc_property_map(sm_selected_map);

  std::map<sm_edge_descriptor,bool> sm_feature_map;
  auto sm_feature_pmap = boost::make_assoc_property_map(sm_feature_map);

  CGAL::IO::read_OM(outname, sm, sm_selected_pmap, sm_feature_pmap);

  std::cout << "vertex selection values:\n";
  for(auto v : vertices(sm)){
    std::cout << std::boolalpha << get(sm_selected_pmap, v) << std::endl;
  }

  std::cout << "edge feature values:\n";
  for(auto e : edges(sm)){
    std::cout  << std::boolalpha << get(sm_feature_pmap, e) << std::endl;
  }
  return 0;
}
