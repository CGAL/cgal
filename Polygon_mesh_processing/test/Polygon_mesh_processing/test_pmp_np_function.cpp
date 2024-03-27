#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <sstream>
#include <iostream>
#include <unordered_map>


namespace CGAL {
template <class PolygonMesh,
          class NamedParameters = parameters::Default_named_parameters >
void my_function_with_named_parameters(PolygonMesh& mesh, const NamedParameters& np = parameters::default_values())
{
  //The class containing all the geometric definitions for the PolygonMesh
  typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type  Traits;
  Traits t;
  CGAL_USE(t);
  //A vertex-index-map that is either taken from the NPs, either an already initialized map for vertex-indices.
  //Also exists for Faces, Edges and Halfedges
  typedef typename CGAL::GetInitializedVertexIndexMap<PolygonMesh, NamedParameters>::type VertexIndexMap;
  //A vertex-point-map either taken from the NPs, either a specified default map.
  typedef typename GetVertexPointMap < PolygonMesh, NamedParameters>::type VPM;

  //The class defining all boost-graph types for the PolygonMesh, like vertex_descriptor and so.
  typedef boost::graph_traits<PolygonMesh> Graph_traits;
  typedef typename Graph_traits::vertex_descriptor vertex_descriptor;

  //in the case no helper function exists, this is how you get a type from a NP
  typedef Static_boolean_property_map<vertex_descriptor, false>                 Default_VCM;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_is_constrained_t,
                                                       NamedParameters,
                                                       Default_VCM>::type       VCM;

  // Normal map
  typedef dynamic_vertex_property_t<typename Traits::Vector_3> Vector_map_tag;
  typedef typename boost::property_map<PolygonMesh, Vector_map_tag>::type Default_vector_map;
  typedef typename internal_np::Lookup_named_param_def<internal_np::vertex_normal_map_t,
                                                       NamedParameters,
                                                       Default_vector_map>::type       VNM;

  using parameters::choose_parameter;
  using parameters::get_parameter;
  using parameters::is_default_parameter;

  //If the NPs provide a vertex-index-map, returns it. Else, returns an initialized vertex-index map.
  VertexIndexMap vim = CGAL::get_initialized_vertex_index_map(mesh, np);
  //If the NPs provide a vertex-point-map, returns it. Else, returns the default boost-graph vpm.
  VPM vpm = choose_parameter(get_parameter(np, internal_np::vertex_point),
                             get_property_map(vertex_point, mesh));
  //boolean NP example. Default value is `false`
  bool do_project = choose_parameter(get_parameter(np, internal_np::do_project), false);

  // If the NPs provide a vertex-normal-map use it, otherwise initialize the default one
  VNM vnm = choose_parameter<Default_vector_map>(get_parameter(np, internal_np::vertex_normal_map),  Vector_map_tag(), mesh);
  if (is_default_parameter<NamedParameters, internal_np::vertex_normal_map_t>::value)
    Polygon_mesh_processing::compute_vertex_normals(mesh, vnm);

  // check is a parameter has been given by the user
  constexpr bool do_project_is_default = is_default_parameter<NamedParameters, internal_np::do_project_t>::value;

  VCM vcm_np = choose_parameter(get_parameter(np, internal_np::vertex_is_constrained), Default_VCM());


  //demonstrates usage for those values.
  for(auto v : vertices(mesh))
  {
    std::cout<<"vertex #"<<get(vim, v)<<" : "<<get(vpm, v)<<" : "<<get(vcm_np, v)<<std::endl;
  }

  if (!do_project_is_default)
  {
    if(do_project)
      std::cout<<"do project"<<std::endl;
    else
      std::cout<<"don't project"<<std::endl;
  }
  else
    std::cout<<"use default for project"<<std::endl;

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

  typedef std::map<boost::graph_traits<SMesh>::vertex_descriptor,bool> VCmap;
  VCmap vcm;
  for(auto v : vertices(sm))
  {
    if ((int)v %2 ==0)
      vcm[v] = true;
    else
      vcm[v] = false;
  }
  typedef boost::associative_property_map<VCmap> Vertex_constrained_pmap;
  Vertex_constrained_pmap vcm_pmap(vcm);
  CGAL::my_function_with_named_parameters(sm);
  CGAL::my_function_with_named_parameters(sm, CGAL::parameters::vertex_is_constrained_map(vcm_pmap)
                                          .do_project(true));
  return 0;
}
