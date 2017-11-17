#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/implicit_shape_smoothing_impl.h>
#include <Eigen/Sparse>

namespace CGAL {


namespace Polygon_mesh_processing {



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
void smooth_modified_curvature_flow(const FaceRange& faces, PolygonMesh& mesh, const NamedParameters& np)
{
  using boost::choose_param;
  using boost::get_param;

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);


  // nb_iterations
  unsigned int nb_iterations = choose_param(get_param(np, internal_np::number_of_iterations), 1);


  std::size_t n = static_cast<int>(vertices(mesh).size());

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

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

  std::cerr << "compute stiffness matrix...";
  smoother.calc_stiff_matrix(stiffness_matrix);
  std::cerr << "done" << std::endl;


  double time = 1e-5;

  for(std::size_t t=0; t<nb_iterations; ++t)
  {
    smoother.setup_system(A, stiffness_matrix, mass_matrix, bx, by, bz, time);
    smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz);
    smoother.update_mesh(Xx, Xy, Xz);
  }

}

// todo: add overloads

// demo API, undocumented
template<typename PolygonMesh>
void setup_mcf_system(PolygonMesh& mesh, Eigen::SparseMatrix<double>& stiffness_matrix)
{

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  std::size_t n = static_cast<int>(vertices(mesh).size());
  stiffness_matrix.resize(n, n);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);
  smoother.calc_stiff_matrix(stiffness_matrix);

}

template<typename PolygonMesh>
void solve_mcf_system(PolygonMesh& mesh, const double& time, unsigned int nb_iter,
	Eigen::SparseMatrix<double>& stiffness_matrix)
{

  // VPmap type
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

  for(unsigned int t = 0; t < nb_iter; ++t)
  {

    smoother.setup_system(A, stiffness_matrix, mass_matrix, bx, by, bz, time);
    smoother.solve_system(A, Xx, Xy, Xz, bx, by, bz);
    smoother.update_mesh(Xx, Xy, Xz);

  }

}






} //Polygon_mesh_processing

} //CGAL
