#include <boost/graph/graph_traits.hpp>
#include <boost/property_map/property_map.hpp>

#include <CGAL/Polygon_mesh_processing/internal/Smoothing/implicit_shape_smoothing_impl.h>

namespace CGAL {


namespace Polygon_mesh_processing {



template<typename PolygonMesh>
void smooth_shape(PolygonMesh& mesh, int nb_iter)
{

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

  Eigen::SparseMatrix<double> stiffness_matrix;
  stiffness_matrix = smoother.calc_stiff_matrix();

  for(unsigned int t=0; t<nb_iter; ++t)
  {
    smoother.solve_system(stiffness_matrix);
  }

}

// demo API, undocumented
template<typename PolygonMesh>
void solve_mcf_system(PolygonMesh& mesh, int nb_iter, Eigen::SparseMatrix<double>& stiffness_matrix)
{

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

  for(unsigned int t=0; t<nb_iter; ++t)
  {
    smoother.solve_system(stiffness_matrix);
  }

}

template<typename PolygonMesh>
void setup_mcf_system(PolygonMesh& mesh, int nb_iter, Eigen::SparseMatrix<double>& stiffness_matrix)
{

  // VPmap type
  typedef typename boost::property_map<PolygonMesh, CGAL::vertex_point_t>::type VertexPointMap;
  VertexPointMap vpmap = get(CGAL::vertex_point, mesh);

  internal::Shape_smoother<PolygonMesh, VertexPointMap> smoother(mesh, vpmap);

  stiffness_matrix = smoother.calc_stiff_matrix();

}




} //Polygon_mesh_processing

} //CGAL
