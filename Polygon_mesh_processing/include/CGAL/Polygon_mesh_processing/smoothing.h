#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H


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




template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void angle_remeshing(PolygonMesh& pmesh, const FaceRange& faces, const NamedParameters& np)
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

    //fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

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
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    //nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    //use weighted angles
    bool use_weights = choose_param(get_param(np, internal_np::use_weights), false);

    internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);
    remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif

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
        remesher.project_to_surface();
    }


#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << "Remeshing done in ";
    std::cout << t.time() << " sec." << std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void angle_remeshing(PolygonMesh& pmesh, const NamedParameters& np)
{
    angle_remeshing(pmesh, faces(pmesh), np);
}

template<typename PolygonMesh>
void angle_remeshing(PolygonMesh& pmesh)
{
    angle_remeshing(pmesh, faces(pmesh), parameters::all_default());
}


template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void area_remeshing(PolygonMesh& pmesh, const FaceRange& faces, const NamedParameters& np)
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

    //fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

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
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    //nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

    //gradient descent precision
    double gd_precision = choose_param(get_param(np, internal_np::gradient_descent_precision), 0.001);

    internal::Compatible_remesher<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            remesher(pmesh, vpmap, vcmap, ecmap);
    remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif

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

        remesher.area_relaxation(gd_precision);
        remesher.project_to_surface();
    }


#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << "Remeshing done in ";
    std::cout << t.time() << " sec." << std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void area_remeshing(PolygonMesh& pmesh, const NamedParameters& np)
{
    area_remeshing(pmesh, faces(pmesh), np);
}

template<typename PolygonMesh>
void area_remeshing(PolygonMesh& pmesh)
{
    area_remeshing(pmesh, faces(pmesh), parameters::all_default());
}

template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void curvature_flow(PolygonMesh& pmesh, const FaceRange& faces, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_REMESHING_VERBOSE
    CGAL::Timer t;
    std::cout << "Remeshing parameters...";
    std::cout.flush();
    t.start();
#endif

    // GeomTraits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    //vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    //fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

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
          internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>
        > ::type ECMap;
    ECMap ecmap = (boost::is_same<ECMap, internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap> >::value)
    ? choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>(pmesh, faces, fimap))
    : choose_param(get_param(np, internal_np::edge_is_constrained),
                   internal::Border_constraint_pmap<PolygonMesh, FaceRange, FIMap>());

    //nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << "\rRemeshing parameters done ("<< t.time() <<" sec)" << std::endl;
    std::cout << "Remesher construction...";
    std::cout.flush();
    t.reset(); t.start();
#endif

    internal::Curvature_flow<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            curvature_remesher(pmesh, vpmap, vcmap, ecmap);
    curvature_remesher.init_remeshing(faces);

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate edges..." << std::endl;
    t.reset(); t.start();
#endif

    curvature_remesher.remove_degenerate_faces();

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

        curvature_remesher.curvature_smoothing();
        //curvature_remesher.good_curvature_smoothing();

        //for now
        //curvature_remesher.project_to_surface();

    }

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << "Remeshing done in ";
    std::cout << t.time() << " sec." << std::endl;
    std::cout<<std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void curvature_flow(PolygonMesh& pmesh, const NamedParameters& np)
{
    curvature_flow(pmesh, faces(pmesh), np);
}

template<typename PolygonMesh>
void curvature_flow(PolygonMesh& pmesh)
{
    curvature_flow(pmesh, parameters::all_default());
}





} // namespace Polygon_mesh_processing
} // namespace CGAL











#endif // CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
