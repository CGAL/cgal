#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

template <typename Traits>
struct My_mesh_1 : public CGAL::Polyhedron_3<Traits> {};

struct My_mesh_2 : public CGAL::Polyhedron_3<Kernel> {};

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
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Polyhedron_3<Traits>
#include <CGAL/boost/graph/graph_traits_inheritance_macros.h>

/// make My_mesh_2 a valid face graph model
// no template parameter, CGAL_GRAPH_TRAITS_INHERITANCE_TEMPLATE_PARAMS is then not defined
#define CGAL_GRAPH_TRAITS_INHERITANCE_CLASS_NAME My_mesh_2
#define CGAL_GRAPH_TRAITS_INHERITANCE_BASE_CLASS_NAME CGAL::Polyhedron_3<Kernel>
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

int main()
{
  typedef My_mesh_1<Kernel>       Mesh1;

  std::vector<Kernel::Point_3> points;
  Mesh1 poly1;
  CGAL::convex_hull_3(points.begin(), points.end(), poly1);

  My_mesh_2 poly2;
  CGAL::convex_hull_3(points.begin(), points.end(), poly2);

  My_mesh_3<Kernel::Point_3> poly3;
  CGAL::convex_hull_3(points.begin(), points.end(), poly3);

  My_mesh_4<3, Kernel> poly4;
  CGAL::convex_hull_3(points.begin(), points.end(), poly4);

  My_mesh_5 poly5;
  CGAL::convex_hull_3(points.begin(), points.end(), poly5);

  return 0;
}
