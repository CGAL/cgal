#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_filtered_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/IO/Color.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_filtered_traits_2<K> GT;

typedef CGAL::Periodic_2_triangulation_vertex_base_2<GT>    Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<CGAL::Color, GT, Vb> VbInfo;

typedef CGAL::Periodic_2_triangulation_face_base_2<GT>      Fb;

typedef CGAL::Triangulation_data_structure_2<VbInfo, Fb>    Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT, Tds>  PDT;

typedef PDT::Point   Point;

int main()
{
  PDT T;

  T.insert(Point(0, 0));
  T.insert(Point(.1, 0));
  T.insert(Point(0, .1));
  T.insert(Point(.2, .2));
  T.insert(Point(.9, 0));

  // Set the color of vertices with degree 6 to red.
  PDT::Vertex_iterator vit;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    if (T.degree(vit) == 6)
      vit->info() = CGAL::RED;

  return 0;
}
