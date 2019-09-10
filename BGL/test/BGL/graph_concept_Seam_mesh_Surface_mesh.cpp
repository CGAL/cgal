#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Simple_cartesian.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor         SM_vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor       SM_halfedge_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor           SM_edge_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor           SM_face_descriptor;

typedef Mesh::Property_map<SM_edge_descriptor, bool>            Seam_edge_pmap;
typedef Mesh::Property_map<SM_vertex_descriptor, bool>          Seam_vertex_pmap;
typedef CGAL::Seam_mesh<Mesh, Seam_edge_pmap, Seam_vertex_pmap> Seam_mesh;

typedef boost::graph_traits< Seam_mesh > Traits;
typedef Traits::edge_descriptor edge_descriptor;
typedef Traits::halfedge_descriptor halfedge_descriptor;
typedef Traits::vertex_descriptor vertex_descriptor;
typedef Traits::face_descriptor face_descriptor;

void concept_check_seam_surface_mesh()
{
  boost::function_requires< boost::GraphConcept<Seam_mesh> >();
  boost::function_requires< boost::VertexListGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::EdgeListGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::IncidenceGraphConcept<Seam_mesh> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Seam_mesh> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Seam_mesh> >();

  boost::function_requires< CGAL::HalfedgeGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::FaceGraphConcept<Seam_mesh> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Seam_mesh> >();

  // null
  boost::graph_traits<Seam_mesh>::null_vertex();
  boost::graph_traits<Seam_mesh>::null_halfedge();
  boost::graph_traits<Seam_mesh>::null_face();
}

int main()
{
  concept_check_seam_surface_mesh();
  return 0;
}
