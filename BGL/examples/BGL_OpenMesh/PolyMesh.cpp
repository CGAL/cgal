#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/property_map.h>
#include <CGAL/Surface_mesh.h>
#include <iostream>


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OMesh;
typedef CGAL::Surface_mesh<Kernel::Point_3> SM;
typedef boost::graph_traits<OMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<OMesh>::halfedge_descriptor halfedge_descriptor;

typedef boost::graph_traits<SM>::vertex_descriptor sm_vertex_descriptor;
typedef boost::graph_traits<SM>::edge_descriptor sm_edge_descriptor;

int main(int argc, char** argv )
{
  OMesh omesh;

  const std::string filename = (argc>1)?argv[1]:CGAL::data_file_path("meshes/cube-meshed.off");

  if (!CGAL::IO::read_polygon_mesh(filename, omesh))
  {
    std::cerr << "Cannot read " << filename << "\n";
    return 1;
  }

  omesh.request_vertex_status();
  omesh.request_edge_status();

  int i = 0;
  for(auto v : vertices(omesh))
  {
    omesh.status(v).set_feature((i%2) == 0);
    ++i;
  }

  i = 0;
  for(auto eit = omesh.edges_begin(); eit != omesh.edges_end(); ++eit)
  {
    omesh.status(*eit).set_feature((i%2) == 0);
    ++i;
  }

  const char* outname= (argc>2)?argv[2]:"out.om";
  OpenMesh::IO::write_mesh(omesh, outname, OpenMesh::IO::Options::Status);

  SM sm;

  std::map<sm_vertex_descriptor,bool> sm_selected_map;
  auto sm_selected_pmap = boost::make_assoc_property_map(sm_selected_map);

  std::map<sm_edge_descriptor,bool> sm_feature_map;
  auto sm_feature_pmap = boost::make_assoc_property_map(sm_feature_map);

  CGAL::IO::read_OM(outname, sm, CGAL::parameters::vertex_is_constrained_map(sm_selected_pmap)
                                                  .edge_is_constrained_map(sm_feature_pmap));


  int nbv=0, nbfv=0;
  for(auto v : vertices(sm)){
    ++nbv;
    if (get(sm_selected_pmap, v)) ++nbfv;
  }
  std::cout << nbfv << " feature vertices out of " <<  nbv << "\n";

  int nbe=0, nbfe=0;
  for(auto e : edges(sm)){
    ++nbe;
    if (get(sm_feature_pmap, e)) ++nbfe;
  }
  std::cout << nbfe << " feature edges out of " <<  nbe << "\n";

  assert(nbv>1 || nbfv!=0);
  assert(nbe>1 || nbfe!=0);

  return 0;
}
