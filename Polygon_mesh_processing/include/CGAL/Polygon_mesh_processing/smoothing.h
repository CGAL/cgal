#ifndef CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H
#define CGAL_POLYGON_MESH_PROCESSING_SMOOTHING_H


#include <CGAL/Polygon_mesh_processing/internal/Smoothing/smoothing_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Smoothing/curvature_flow_impl.h>
#include <CGAL/Polygon_mesh_processing/internal/Isotropic_remeshing/remesh_impl.h>

#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

#ifdef CGAL_PMP_REMESHING_VERBOSE
#include <CGAL/Timer.h>
#endif

namespace CGAL {

namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* @brief smoothes a triangulated region of a polygon mesh.
* This function imrpoves the angles between triangle edges by moving not
* constrained vertices so that each pair of adjacent angles becomes equal.
* Optionally, small angles may carry more weight than larger ones. Projection
* to the initial surface is performed as a last step.
* for better result.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be remeshed
* @param faces the range of triangular faces defining one or several surface patches to be remeshed
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of atomic operations performed (listed in the above description)
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. A constrained edge can be split
*    or collapsed, but not flipped, nor its endpoints moved by smoothing.
*    Note that patch boundary edges (i.e. incident to only one face in the range)
*    are always considered as constrained edges.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during remeshing.
*  \cgalParamEnd
*  \cgalParamBegin{use_weights} If `true`, small angles carry more weight than larger ones.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
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
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.init_remeshing(faces);


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

/*!
* \ingroup PMP_meshing_grp
* @brief uses triangle area as a criterion for smoothing.
* This function imrpoves the overall distribution of points over a mesh area
* by trying to form triangles as uniform as possible. Should be used in combination with angle remeshing
* to avoid creation of long and skiny triangles. Vertices are moved towards equalizing adjacent triangle
* areas using gradient descent.
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be remeshed
* @param faces the range of triangular faces defining one or several surface patches to be remeshed
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of atomic operations performed (listed in the above description)
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. A constrained edge can be split
*    or collapsed, but not flipped, nor its endpoints moved by smoothing.
*    Note that patch boundary edges (i.e. incident to only one face in the range)
*    are always considered as constrained edges.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during remeshing.
*  \cgalParamEnd
*  \cgalParamBegin{gradient_descent_precision} The precision which is met during gradient descent refers to
*    the relative energy between iterations of each triangle element which is minimized
*    while one of its vertices is being moved. Triangle energy is defined based on its area compared to
*    the average area of all triangles adjacent to the vertex that is being moved.  Defaults to 0.001.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
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
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    remesher.init_remeshing(faces);

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

/*!
* \ingroup PMP_meshing_grp
* @brief todo
*
* @tparam PolygonMesh model of `MutableFaceGraph`.
*         The descriptor types `boost::graph_traits<PolygonMesh>::%face_descriptor`
*         and `boost::graph_traits<PolygonMesh>::%halfedge_descriptor` must be
*         models of `Hashable`.
*         If `PolygonMesh` has an internal property map for `CGAL::face_index_t`,
*         and no `face_index_map` is given
*         as a named parameter, then the internal one should be initialized
* @tparam FaceRange range of `boost::graph_traits<PolygonMesh>::%face_descriptor`,
          model of `Range`. Its iterator type is `ForwardIterator`.
* @tparam NamedParameters a sequence of \ref namedparameters
*
* @param pmesh a polygon mesh with triangulated surface patches to be remeshed
* @param faces the range of triangular faces defining one or several surface patches to be remeshed
* @param np optional sequence of \ref namedparameters among the ones listed below
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of atomic operations performed (listed in the above description)
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. A constrained edge can be split
*    or collapsed, but not flipped, nor its endpoints moved by smoothing.
*    Note that patch boundary edges (i.e. incident to only one face in the range)
*    are always considered as constrained edges.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during remeshing.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
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
#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif
    curvature_remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_REMESHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif
    curvature_remesher.init_remeshing(faces);

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

        curvature_remesher.curvature_smoothing(); // normalized version sec 5.5
        //curvature_remesher.good_curvature_smoothing(); // formula 14

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
