#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H


#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/smoothing_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>




namespace CGAL {

namespace Polygon_mesh_processing {





template<typename PolygonMesh, typename NamedParameters, typename FaceRange, typename EdgeRange>
void angle_remeshing(PolygonMesh& pmesh, const FaceRange& faces, const EdgeRange& edges, const NamedParameters& np)
{

    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;



    // extract some edges
    /*
    typedef typename boost::graph_traits<PolygonMesh>::edge_iterator edge_iterator;
    typedef std::pair<edge_iterator,edge_iterator> p_edges;

    typename boost::graph_traits<PolygonMesh>::edge_iterator ei, ei_end;


    for(boost::tie(ei, ei_end) = edges; ei != ei_end; ++ei)
    {
        //std::cout<<"p: "<<p<<std::endl;
        std::cout<<source(*ei, pmesh)<<"->"<<target(*ei, pmesh)<<std::endl;
    }
    */


    //vcmap - not sure yet
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());



    //create ecmap with Edge_constraint_map struct
    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

    //ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Edge_constraint_map<PolygonMesh, EdgeRange>
        > ::type ECMap;
    ECMap ecmap = choose_param(get_param(np, internal_np::edge_is_constrained),
                               internal::Edge_constraint_map<PolygonMesh, EdgeRange>()); // pass constrained edges range in this constructor


    CGAL::Polygon_mesh_processing::internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits> remesher(pmesh, vpmap, vcmap, ecmap);
    remesher.init_remeshing(faces);
    remesher.angle_relaxation();
    remesher.project_to_surface();


}



template<typename PolygonMesh, typename NamedParameters, typename FaceRange, typename EdgeRange>
void area_remeshing(PolygonMesh& pmesh,  const FaceRange& faces, const EdgeRange& edges, const NamedParameters& np)
{

    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;


    //bool do_project = choose_param(get_param(np, internal_np::do_project), true);



    //vcmap - not sure yet
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());



    //create ecmap with Edge_constraint_map struct
    typedef typename boost::graph_traits<PolygonMesh>::edge_descriptor edge_descriptor;

    //ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Edge_constraint_map<PolygonMesh, EdgeRange>
        > ::type ECMap;
    ECMap ecmap = choose_param(get_param(np, internal_np::edge_is_constrained),
                               internal::Edge_constraint_map<PolygonMesh, EdgeRange>());



    CGAL::Polygon_mesh_processing::internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits> remesher(pmesh, vpmap, vcmap, ecmap);
    remesher.init_remeshing(faces);
    remesher.area_relaxation();


}










} // namespace Polygon_mesh_processing
} // namespace CGAL











#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
