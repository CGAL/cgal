#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_filtered_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/IO/Color.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_filtered_traits_3<K> GT;

typedef CGAL::Periodic_3_triangulation_ds_vertex_base_3<> VbDS;
typedef CGAL::Triangulation_vertex_base_3<GT,VbDS> Vb;

typedef CGAL::Periodic_3_triangulation_ds_cell_base_3<> CbDS;
typedef CGAL::Triangulation_cell_base_3<GT,CbDS> Cb;

typedef CGAL::Triangulation_vertex_base_with_info_3<CGAL::Color, GT, Vb> VbInfo;
typedef CGAL::Triangulation_data_structure_3<VbInfo, Cb>    TDS;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT, TDS>  PDT;

typedef PDT::Point   Point;

int main()
{
  PDT T;

  T.insert(Point(0,0,0));
  T.insert(Point(.1,0,0));
  T.insert(Point(0,.1,0));
  T.insert(Point(0,0,.1));
  T.insert(Point(.2,.2,.2));
  T.insert(Point(.9,0,.1));

  // Set the color of finite vertices of degree 6 to red.
  PDT::Vertex_iterator vit;
  for (vit = T.vertices_begin(); vit != T.vertices_end(); ++vit)
    if (T.degree(vit) == 6)
      vit->info() = CGAL::RED;

  return 0;
}
