#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H


#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/smoothing_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>




namespace CGAL {

namespace Polygon_mesh_processing {





template<typename PolygonMesh, typename NamedParameters, typename FaceRange>
void angle_remeshing(PolygonMesh& pmesh, const FaceRange& faces, const NamedParameters& np)
{

    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type Traits;


    CGAL::Polygon_mesh_processing::internal::Angle_remesher<PolygonMesh, VertexPointMap, Traits> remesher(pmesh, vpmap);
    remesher.init_remeshing(faces);
    remesher.angle_relaxation();
    remesher.project_to_surface();


}



template<typename PolygonMesh, typename NamedParameters, typename FaceRange>
void area_remeshing(PolygonMesh& pmesh,  const FaceRange& faces, const NamedParameters& np)
{

    using boost::choose_param;
    using boost::get_param;

    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::const_type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_const_property_map(CGAL::vertex_point, pmesh));


    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;


    //bool do_project = choose_param(get_param(np, internal_np::do_project), true);


    CGAL::Polygon_mesh_processing::internal::Area_remesher<PolygonMesh, VertexPointMap, GeomTraits> remesher(pmesh, vpmap);
    remesher.init_remeshing(faces);
    remesher.area_relaxation();


}










} // namespace Polygon_mesh_processing
} // namespace CGAL











#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
