#include <CGAL/Polygon_mesh_processing/internal/named_function_params.h>
#include <CGAL/Polygon_mesh_processing/internal/named_params_helper.h>

namespace CGAL {
namespace Polygon_mesh_processing {
namespace internal {

template <typename PolygonMesh, typename Descriptor, typename Triangle>
void construct_triangle(const Descriptor& f, const PolygonMesh& mesh_, Triangle& triangle)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<PolygonMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type Vpm;
  Vpm vpmap_ = get(boost::vertex_point, mesh_);

  halfedge_descriptor h = halfedge(f, mesh_);
  vertex_descriptor v1 = target(h, mesh_);
  vertex_descriptor v2 = target(next(h, mesh_), mesh_);
  vertex_descriptor v3 = target(next(next(h, mesh_), mesh_), mesh_);
  triangle = Triangle(get(vpmap_, v1), get(vpmap_, v2), get(vpmap_, v3));
}

template <typename PolygonMesh, typename Descriptor>
double sqlength(const Descriptor& v1, const Descriptor& v2, const PolygonMesh& mesh_)
{
  typedef typename boost::property_map<PolygonMesh, boost::vertex_point_t>::const_type Vpm;
  Vpm vpmap_ = get(boost::vertex_point, mesh_);
  return to_double(CGAL::squared_distance(get(vpmap_, v1), get(vpmap_, v2)));
}

template <typename PolygonMesh, typename Descriptor>
double sqlength(const Descriptor& h, const PolygonMesh& mesh_)
{
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  vertex_descriptor v1 = target(h, mesh_);
  vertex_descriptor v2 = source(h, mesh_);
  return sqlength(v1, v2, mesh_);
}

}}}
