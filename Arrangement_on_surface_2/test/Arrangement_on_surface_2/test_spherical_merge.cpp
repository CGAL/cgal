// #define CGAL_IDENTIFICATION_XY  2

#include <string>
#include <vector>

#include <CGAL/Gmpq.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/IO/Arr_iostream.h>

typedef CGAL::Gmpq                                           Number_type;
typedef CGAL::Cartesian<Number_type>                         Kernel;
typedef CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>    Geom_traits_2;
typedef Geom_traits_2::Point_2                               Point_2;
typedef Geom_traits_2::Curve_2                               Curve_2;
typedef Geom_traits_2::X_monotone_curve_2                    X_monotone_curve_2;
typedef CGAL::Arr_spherical_topology_traits_2<Geom_traits_2> Topol_traits_2;
typedef CGAL::Arrangement_on_surface_2<Geom_traits_2, Topol_traits_2>
                                                             Arrangement_2;

typedef Arrangement_2::Halfedge_handle                       Halfedge_handle;
typedef Arrangement_2::Vertex_handle                         Vertex_handle;

int main()
{
  Arrangement_2 arr;

  Point_2 p1(0, 0, -1), p2(1, 0, 0), p3(0, 0, 1), p4(0, -1, 0);
  X_monotone_curve_2 c1 = X_monotone_curve_2(p1, p2);
  X_monotone_curve_2 c2 = X_monotone_curve_2(p2, p3);
  X_monotone_curve_2 c3 = X_monotone_curve_2(p3, p4);
  X_monotone_curve_2 c4 = X_monotone_curve_2(p4, p1);
  Halfedge_handle he1 = arr.insert_in_face_interior(c1, arr.reference_face());
  Vertex_handle v1 = he1->source();
  Vertex_handle v2 = he1->target();
  //std::cout << v1->point() << std::endl;
  //std::cout << v2->point() << std::endl;
  Halfedge_handle he2 = arr.insert_from_left_vertex(c2, v2);
  Vertex_handle v3 = he2->target();
  //std::cout << v3->point() << std::endl;
  Halfedge_handle he3 = arr.insert_from_right_vertex(c3, v3);
  Vertex_handle v4 = he3->target();
  //std::cout << v4->point() << std::endl;
  Halfedge_handle e4 = arr.insert_at_vertices(c4, v4, v1);
  //std::cout << arr << std::endl;

  CGAL_assertion(v1->degree() == 2);
  CGAL_assertion(v2->degree() == 2);
  CGAL_assertion(v3->degree() == 2);
  CGAL_assertion(v4->degree() == 2);

  const Geom_traits_2* traits = arr.geometry_traits();

  if (!traits->are_mergeable_2_object()(he1->curve(), he1->next()->curve()))
    return -1;
  X_monotone_curve_2 c12;
  traits->merge_2_object()(he1->curve(), he1->next()->curve(), c12);
  arr.merge_edge(he1, he1->next(), c12);

  if (!traits->are_mergeable_2_object()(he3->curve(), he3->next()->curve()))
    return -1;
  X_monotone_curve_2 c34;
  traits->merge_2_object()(he3->curve(), he3->next()->curve(), c34);
  arr.merge_edge(he3, he3->next(), c34);

  std::cout << arr << std::endl;
  return 0;
}
