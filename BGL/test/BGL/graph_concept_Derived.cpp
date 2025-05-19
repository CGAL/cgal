#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>

#include <boost/graph/graph_concepts.hpp>
#include <CGAL/boost/graph/graph_concepts.h>


typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;

template <typename Traits>
struct My_mesh_1 : public CGAL::Polyhedron_3<Traits, CGAL::Polyhedron_items_with_id_3> {};

struct My_mesh_2 : public CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> {};

template <typename PT>
struct My_mesh_3 : public CGAL::Surface_mesh<PT> {};

struct My_mesh_5 : public CGAL::Surface_mesh<Kernel::Point_3> {};

// dim could be hard-coded but for the purpose of the example it is left
template <int dim, typename K>
struct My_mesh_4 :
  CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
         <2, dim, CGAL::Linear_cell_complex_traits<dim, K> >::type
{};

/// make My_mesh_1 a valid face graph model
#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS typename Traits
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_1<Traits>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Polyhedron_3<Traits, CGAL::Polyhedron_items_with_id_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

/// make My_mesh_2 a valid face graph model
// no template parameter, CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS is then not defined
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_2
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

/// make My_mesh_3 a valid face graph model
#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS typename PT
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_3<PT>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<PT>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

/// make My_mesh_4 a valid face graph model
#define CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS int dim, typename K
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_4<dim, K>
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME typename CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper\
         <2, dim, CGAL::Linear_cell_complex_traits<dim, K> >::type
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

/// make My_mesh_5 a valid face graph model
// no template parameter, CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS is then not defined
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_5
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Surface_mesh<Kernel::Point_3>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>


template <class Graph>
void concept_check()
{
  typedef boost::graph_traits< Graph > Traits;
  typedef typename Traits::edge_descriptor edge_descriptor;
  typedef typename Traits::halfedge_descriptor halfedge_descriptor;
  typedef typename Traits::vertex_descriptor vertex_descriptor;
  typedef typename Traits::face_descriptor face_descriptor;

  boost::function_requires< boost::GraphConcept<Graph> >();
  boost::function_requires< boost::VertexListGraphConcept<Graph> >();
  boost::function_requires< boost::EdgeListGraphConcept<Graph> >();
  boost::function_requires< boost::IncidenceGraphConcept<Graph> >();
  boost::function_requires< boost::AdjacencyMatrixConcept<Graph> >();
  boost::function_requires< boost::BidirectionalGraphConcept<Graph> >();
  boost::function_requires< CGAL::HalfedgeGraphConcept<Graph> >();
  boost::function_requires< CGAL::HalfedgeListGraphConcept<Graph> >();
  boost::function_requires< CGAL::FaceGraphConcept<Graph> >();
  boost::function_requires< CGAL::FaceListGraphConcept<Graph> >();
  boost::function_requires< CGAL::MutableHalfedgeGraphConcept<Graph> >();
  boost::function_requires< CGAL::MutableFaceGraphConcept<Graph> >();

  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, halfedge_descriptor, CGAL::halfedge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, edge_descriptor, boost::edge_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, edge_descriptor, boost::edge_weight_t> >();
  boost::function_requires< boost::concepts::PropertyGraph<
    Graph, vertex_descriptor, CGAL::vertex_point_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, vertex_descriptor, boost::vertex_index_t> >();
  boost::function_requires< boost::concepts::ReadablePropertyGraph<
    Graph, face_descriptor, CGAL::face_index_t> >();

  // null
  boost::graph_traits<Graph>::null_vertex();
  boost::graph_traits<Graph>::null_halfedge();
  boost::graph_traits<Graph>::null_face();
}

int main()
{
  concept_check<My_mesh_1<Kernel>>();
  concept_check<My_mesh_2>();
  concept_check<My_mesh_3<Point_3> >();
  concept_check<My_mesh_4<3,Kernel>>();
  concept_check<My_mesh_5>();
  return 0;
}
