#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/implicit_shape_smoothing_impl.h>
#include <Eigen/Sparse>

namespace CGAL {
namespace Polygon_mesh_processing {

/*!
* \ingroup PMP_meshing_grp
* smooths the overall shape of the mesh by using the mean curvature flow.
* The effect depends only on the curvature of each area.
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
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref namedparameters among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Exact constructions kernels are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations for the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. Vertices that belong to constrained
*    edges are not modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_curvature_flow(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
    using boost::choose_param;
    using boost::get_param;

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    CGAL::Timer t;
    std::cout << "Smoothing parameters...";
    std::cout.flush();
    t.start();
#endif

    // GeomTraits
    typedef typename GetGeomTraits<PolygonMesh, NamedParameters>::type GeomTraits;

    // vpmap
    typedef typename GetVertexPointMap<PolygonMesh, NamedParameters>::type VertexPointMap;
    VertexPointMap vpmap = choose_param(get_param(np, internal_np::vertex_point),
                                 get_property_map(CGAL::vertex_point, pmesh));

    // fimap
    typedef typename GetFaceIndexMap<PolygonMesh, NamedParameters>::type FIMap;
    FIMap fimap = choose_param(get_param(np, internal_np::face_index),
                             get_property_map(face_index, pmesh));

    // vcmap
    typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
    typedef typename boost::lookup_named_param_def <
        internal_np::vertex_is_constrained_t,
        NamedParameters,
        internal::No_constraint_pmap<vertex_descriptor>//default
      > ::type VCMap;
    VCMap vcmap = choose_param(get_param(np, internal_np::vertex_is_constrained),
                               internal::No_constraint_pmap<vertex_descriptor>());

    // ecmap
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

    // nb_iterations
    unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << "\rSmoothing parameters done ("<< t.time() <<" sec)" << std::endl;
    std::cout << "Remesher construction...";
    std::cout.flush();
    t.reset(); t.start();
#endif

    internal::Curvature_flow<PolygonMesh, VertexPointMap, VCMap, ECMap, GeomTraits>
            curvature_remesher(pmesh, vpmap, vcmap, ecmap);
#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Removing degenerate faces..." << std::endl;
    t.reset(); t.start();
#endif

    curvature_remesher.remove_degenerate_faces();

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "Initializing..." << std::endl;
    t.reset(); t.start();
#endif

    curvature_remesher.init_smoothing(faces);

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << " done ("<< t.time() <<" sec)." << std::endl;
    std::cout << "#iter = " << nb_iterations << std::endl;
    std::cout << "Shape smoothing..." << std::endl;
    t.reset(); t.start();
#endif

    for(unsigned int i=0; i<nb_iterations; ++i)
    {

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
        std::cout << " * Iteration " << (i + 1) << " *" << std::endl;
#endif

    curvature_remesher.curvature_smoothing();

    }

#ifdef CGAL_PMP_SMOOTHING_VERBOSE
    t.stop();
    std::cout << "Shape smoothing done in ";
    std::cout << t.time() << " sec." << std::endl;
    std::cout<<std::endl;
#endif

}

template<typename PolygonMesh, typename NamedParameters>
void smooth_curvature_flow(PolygonMesh& pmesh, const NamedParameters& np)
{
    smooth_curvature_flow(faces(pmesh), pmesh, np);
}

template<typename PolygonMesh>
void smooth_curvature_flow(PolygonMesh& pmesh)
{
    smooth_curvature_flow(pmesh, parameters::all_default());
}

/*!
* \ingroup PMP_meshing_grp
* smooths the overall shape of the mesh by using a modified algorithm based on the mean curvature flow.
* The effect depends only on the curvature of each area and the convergence is achieved to a conformal map of the initial surface
* to a sphere.
*
* @todo add time parameter.
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
* @param pmesh a polygon mesh with triangulated surface patches to be smoothed.
* @param faces the range of triangular faces defining one or several surface patches to be smoothed.
* @param np optional sequence of \ref namedparameters among the ones listed below.
*
* \cgalNamedParamsBegin
*  \cgalParamBegin{geom_traits} a geometric traits class instance, model of `Kernel`.
*    Kernels with exact constructions are not supported by this function.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_point_map} the property map with the points associated
*    to the vertices of `pmesh`. Instance of a class model of `ReadWritePropertyMap`.
*  \cgalParamEnd
*  \cgalParamBegin{number_of_iterations} the number of iterations of the
*    sequence of the smoothing iterations performed.
*  \cgalParamEnd
*  \cgalParamBegin{face_index_map} a property map containing the index of each face of `pmesh`.
*  \cgalParamEnd
*  \cgalParamBegin{edge_is_constrained_map} a property map containing the
*    constrained-or-not status of each edge of `pmesh`. Vertices that belong to constrained
*    edges are not modified at all during smoothing.
*  \cgalParamEnd
*  \cgalParamBegin{vertex_is_constrained_map} a property map containing the
*    constrained-or-not status of each vertex of `pmesh`. A constrained vertex
*    cannot be modified at all during smoothing.
*  \cgalParamEnd
* \cgalNamedParamsEnd
*/
template<typename PolygonMesh, typename FaceRange, typename NamedParameters>
void smooth_modified_curvature_flow(const FaceRange& faces, PolygonMesh& pmesh, const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, pmesh);


  // nb_iterations
  unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);


  std::size_t n = static_cast<int>(vertices(pmesh).size());

  typedef typename Eigen::VectorXd Eigen_vector;
  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;

  Eigen_matrix A(n, n);
  Eigen_matrix stiffness_matrix(n, n);
  Eigen_matrix mass_matrix(n, n);
  Eigen_vector bx(n);
  Eigen_vector by(n);
  Eigen_vector bz(n);
  Eigen_vector Xx(n);
  Eigen_vector Xy(n);
  Eigen_vector Xz(n);

  // use resize

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(pmesh, vpmap);




  double time = 1e-5;

  for(std::size_t t=0; t<nb_iterations; ++t)
  {

    std::cerr << "compute stiffness matrix...";
    smoother.calc_stiff_matrix(stiffness_matrix);
    // print out stiff matrix
    std::ofstream out("data/stiff_matrix.dat");
    out<<stiffness_matrix<<std::endl;
    out.close();


    smoother.setup_system(A, stiffness_matrix, mass_matrix, bx, by, bz, time);

    /*std::ofstream out("data/mass_matrix.dat");
    out<<mass_matrix<<std::endl;
    out.close();*/

    smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz);
    smoother.update_mesh(Xx, Xy, Xz);
  }

}

// todo: add overloads

// demo API, undocumented
template<typename PolygonMesh>
void setup_mcf_system(PolygonMesh& mesh, Eigen::SparseMatrix<double>& stiffness_matrix)
{
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  std::size_t n = static_cast<int>(vertices(mesh).size());
  stiffness_matrix.resize(n, n);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);
  smoother.calc_stiff_matrix(stiffness_matrix);
}

template<typename PolygonMesh>
void solve_mcf_system(PolygonMesh& mesh, const double& time,
                      Eigen::SparseMatrix<double>& stiffness_matrix)
{
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  typedef typename Eigen::VectorXd Eigen_vector;
  typedef typename Eigen::SparseMatrix<double> Eigen_matrix;

  std::size_t n = static_cast<int>(vertices(mesh).size());

  // temp
  Eigen_matrix A(n, n);
  Eigen_matrix mass_matrix(n, n);
  Eigen_vector bx(n);
  Eigen_vector by(n);
  Eigen_vector bz(n);
  Eigen_vector Xx(n);
  Eigen_vector Xy(n);
  Eigen_vector Xz(n);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

  smoother.setup_system(A, stiffness_matrix, mass_matrix, bx, by, bz, time);
  smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz);
  smoother.update_mesh(Xx, Xy, Xz);
}


} //Polygon_mesh_processing
} //CGAL
