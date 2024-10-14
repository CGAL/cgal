#pragma once

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/property_map.h>
#include <iostream>
#include <fstream>
#include <map>

namespace CGAL {
namespace IO {

template <typename SM, typename VSelectionPM, typename EFeaturePM>
bool read_OM(std::string fname, SM& sm, VSelectionPM vspm, EFeaturePM efpm)
{
  typedef OpenMesh::PolyMesh_ArrayKernelT<> OMesh;
  typedef boost::graph_traits<OMesh>::vertex_descriptor om_vertex_descriptor;
  typedef boost::graph_traits<SM>::vertex_descriptor sm_vertex_descriptor;
  typedef boost::graph_traits<OMesh>::halfedge_descriptor om_halfedge_descriptor;
  typedef boost::graph_traits<SM>::halfedge_descriptor sm_halfedge_descriptor;

  OMesh omesh;
  OpenMesh::IO::Options options = OpenMesh::IO::Options::Status;
  bool ok = OpenMesh::IO::read_mesh(omesh, fname, options);
  if(! ok){
    return false;
  }

  std::map<om_vertex_descriptor,sm_vertex_descriptor> v2v;
  auto v2vpmap = boost::make_assoc_property_map(v2v);

  std::map<om_halfedge_descriptor,sm_halfedge_descriptor> h2h;
  auto h2hpmap = boost::make_assoc_property_map(h2h);

  CGAL::copy_face_graph<OMesh,SM>(omesh, sm, CGAL::parameters::vertex_to_vertex_map(v2vpmap).halfedge_to_halfedge_map(h2hpmap));

  if(options.vertex_has_status()){
    for(auto v : vertices(omesh)){
        put(vspm, v2v[v], omesh.status(v).selected());
    }
  }else{
    std::cout << "no vertex status" << std::endl;
  }

  if(options.edge_has_status()){
    for(auto e : edges(omesh)){
        auto sme = edge(h2h[halfedge(e,omesh)], sm);
        put(efpm, sme , omesh.status(OpenMesh::EdgeHandle(e.idx())).feature());
    }
  }else{
    std::cout << "no edge status" << std::endl;
  }
  return true;
}

} // namespace IO
} // namespace CGAL


