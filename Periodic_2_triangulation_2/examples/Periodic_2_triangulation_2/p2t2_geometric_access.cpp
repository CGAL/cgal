#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_traits_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_triangulation_traits_2<K>          PK;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<PK>       P2DT2;

typedef PK::Point_2        Point;
typedef PK::Triangle_2     Triangle;

typedef P2DT2::Periodic_triangle           Periodic_triangle;
typedef P2DT2::Periodic_triangle_iterator  Periodic_triangle_iterator;
typedef P2DT2::Iterator_type               Iterator_type;

int main()
{
  P2DT2 T;

  T.insert(Point(0, 0));
  T.insert(Point(0, 0.5));
  T.insert(Point(0.5, 0));

  Periodic_triangle pt;
  Triangle t_bd;

  // Extracting the triangles that have a non-empty intersection with
  // the original domain of the 1-sheeted covering space
  for (Periodic_triangle_iterator ptit = T.periodic_triangles_begin(P2DT2::UNIQUE_COVER_DOMAIN);
       ptit != T.periodic_triangles_end(P2DT2::UNIQUE_COVER_DOMAIN); ++ptit)
    {
      pt = *ptit;
      if (! (pt[0].second.is_null() && pt[1].second.is_null() && pt[2].second.is_null()) )
        {
          // Convert the current Periodic_triangle to a Triangle if it is
          // not strictly contained inside the original domain.
          // Note that this requires EXACT constructions to be exact!
          t_bd = T.triangle(pt);
        }
    }
}
