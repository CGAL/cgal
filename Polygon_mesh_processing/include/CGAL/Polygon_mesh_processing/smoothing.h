#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H

#define CGAL_PMP_REMESHING_VERBOSE


#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/smoothing_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/curvature_flow_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef CGAL_PMP_REMESHING_VERBOSE
#include <CGAL/Timer.h>
#endif





namespace CGAL {

namespace Polygon_mesh_processing {





template<typename PolygonMesh, typename NamedParameters, typename FaceRange, typename EdgeRange>
void compatible_remeshing(PolygonMesh& pmesh, const FaceRange& faces, const EdgeRange& edges, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_REMESHING_VERBOSE
  CGAL::Timer t;
  std::cout << "Remeshing parameters...";
  std::cout.flush();
  t.start();
#endif

    //geom traits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    //vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    //vcmap
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());

    //ecmap
    typedef typename boost::lookup_named_param_def <
          internal_np::edge_is_constrained_t,
          NamedParameters,
          internal::Edge_constraint_map<PolygonMesh, EdgeRange>
        > ::type ECMap;
    // either fill map with given constrined edges or use default constructor for an NULL map
    ECMap ecmap = (boost::is_same<ECMap, internal::Edge_constraint_map<PolygonMesh, EdgeRange> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Edge_constraint_map<PolygonMesh, EdgeRange>(edges))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Edge_constraint_map<PolygonMesh, EdgeRange>());

    //nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    //gradient descent precision
    double gd_precision = choose_param(get_param(np, internal_np::gradient_descent_precision), 0.001);

    //use weighted angles
    bool use_weights = choose_param(get_param(np, internal_np::use_weights), false);


#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "\rRemeshing parameters done ("<< t.time() <<" sec)" << std::endl;
  std::cout << "Remesher construction...";
  std::cout.flush();
  t.reset(); t.start();
#endif

    CGAL::Polygon_mesh_processing::internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);
    remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Removing degenerate edges and faces..." << std::endl;
  t.reset(); t.start();
#endif

    //remesher.collapse_short_edges();
    remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "#iter = " << nb_iterations << std::endl;
  std::cout << "Remeshing ..." << std::endl;
  t.reset(); t.start();
#endif

    for(unsigned int i=0; i<nb_iterations; ++i)
    {
#ifdef CGAL_PMP_REMESHING_VERBOSE
    std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif
        remesher.angle_relaxation(use_weights);
        remesher.area_relaxation(gd_precision);
        remesher.angle_relaxation(use_weights);
        remesher.project_to_surface();
    }


#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "Remeshing done in ";
  std::cout << t.time() << " sec." << std::endl;
  std::cout<<std::endl;
#endif

}

template<typename PolygonMesh>
void compatible_remeshing(PolygonMesh& pmesh)
{
    compatible_remeshing(pmesh, faces(pmesh), edges(pmesh), parameters::all_default());
}

template<typename PolygonMesh, typename NamedParameters>
void compatible_remeshing(PolygonMesh& pmesh, const NamedParameters& np)
{
    compatible_remeshing(pmesh, faces(pmesh), edges(pmesh), np);
}



template<typename PolygonMesh, typename NamedParameters>
void curvature_flow(PolygonMesh& pmesh, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_REMESHING_VERBOSE
  CGAL::Timer t;
  std::cout << "Remeshing parameters...";
  std::cout.flush();
  t.start();
#endif

    //vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    // GeomTraits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;


#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "\rRemeshing parameters done ("<< t.time() <<" sec)" << std::endl;
  std::cout << "Remesher construction...";
  std::cout.flush();
  t.reset(); t.start();
#endif

    internal::Curvature_flow<PolygonMesh, VertexPointMap, GeomTraits> curvature_remesher(pmesh, vpmap);
    curvature_remesher.init_remeshing(faces(pmesh));

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Removing degenerate edges..." << std::endl;
  t.reset(); t.start();
#endif

    //curvature_remesher.collapse_short_edges();
    curvature_remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << " done ("<< t.time() <<" sec)." << std::endl;
  std::cout << "Remeshing ..." << std::endl;
  t.reset(); t.start();
#endif

    curvature_remesher.curvature_smoothing();
    curvature_remesher.project_to_surface();


#ifdef CGAL_PMP_REMESHING_VERBOSE
  t.stop();
  std::cout << "Remeshing done in ";
  std::cout << t.time() << " sec." << std::endl;
  std::cout<<std::endl;
#endif


}

template<typename PolygonMesh>
void curvature_flow(PolygonMesh& pmesh)
{
    curvature_flow(pmesh, parameters::all_default());
}




} // namespace Polygon_mesh_processing
} // namespace CGAL











#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
