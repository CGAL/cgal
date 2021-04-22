#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <sstream>
#include <iostream>
#include <unordered_map>

#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/boost/graph/Named_function_parameters.h>

namespace CGAL {
template <class PolygonMesh,
          class NamedParameters>
void np_function(PolygonMesh& mesh, const NamedParameters& np)
{

  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type  Traits;
  typedef typename CGAL::GetInitializedVertexIndexMap<PolygonMesh, NamedParameters>::type VertexIndexMap;
  typedef typename GetVertexPointMap < PolygonMesh, NamedParameters>::type VPM;
  typedef typename GetFaceNormalMap < PolygonMesh, NamedParameters>::type FNM;

  typedef typename Traits::Point_3 Point_3;
  typedef typename Traits::Vector_3 Vector_3;

  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;
  typedef typename Graph_traits::edge_descriptor edge_descriptor;
  typedef typename Graph_traits::halfedge_descriptor halfedge_descriptor;
  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;


  FNM default_fvmap;
  VertexIndexMap vim = CGAL::get_initialized_vertex_index_map(mesh, np);
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, mesh));
  FNM fnm = choose_parameter(get_parameter(np, internal_np::face_normal),
                             default_fvmap);

  bool has_fnm = !is_default_parameter(get_parameter(np, internal_np::face_normal));

  if(has_fnm)
    for(auto f : faces(mesh))
    {
      std::cout<<get(fnm, f)<<std::endl;
    }
  for(auto v : vertices(mesh))
  {
    std::cout<<"vertex #"<<get(vim, v)<<" : "<<get(vpm, v)<<std::endl;
  }

}
}

int main()
{
  typedef CGAL::Surface_mesh<CGAL::Epick::Point_3> SMesh;
  SMesh sm;
  std::stringstream poly("OFF\n"
                         "4 4 0\n"
                         "0 0 0\n"
                         "0 0 1\n"
                         "0 1 0\n"
                         "1 0 0\n"
                         "3  3 1 2\n"
                         "3  0 1 3\n"
                         "3  0 3 2\n"
                         "3  0 2 1\n");
  poly >> sm;
  typedef std::map<boost::graph_traits<SMesh>::face_descriptor,CGAL::Epick::Vector_3> FNmap;
  FNmap fnm;
  for(auto f : faces(sm))
    fnm[f] = {0,0,1};
  typedef boost::associative_property_map<FNmap> Face_normal_pmap;
  Face_normal_pmap fn_pmap(fnm);
  CGAL::np_function(sm, CGAL::parameters::all_default());
  CGAL::np_function(sm, CGAL::parameters::face_normal_map(fn_pmap));
  return 0;
}
