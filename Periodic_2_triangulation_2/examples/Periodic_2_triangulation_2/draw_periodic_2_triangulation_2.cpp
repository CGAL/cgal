#include <CGAL/Periodic_2_Delaunay_triangulation_2.h>
#include <CGAL/Periodic_2_Delaunay_triangulation_traits_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/draw_periodic_2_triangulation_2.h>

#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_2_Delaunay_triangulation_traits_2<K> GT;
typedef CGAL::Periodic_2_Delaunay_triangulation_2<GT>       PDT;

typedef PDT::Point                                          Point;

int main(int argc, char* argv[])
{
  // Declare periodic triangulation 2D
  PDT T;

  // Read points and insert in T
  Point p;
  std::ifstream ifs((argc > 1) ? argv[1] : "data/data1.dt.cin");
  if (ifs)
  {
    while (ifs >> p)
    { T.insert(p); }
    CGAL_assertion(T.is_valid());

    if( T.is_triangulation_in_1_sheet())
    { T.convert_to_9_sheeted_covering(); }

    // Draw the periodic triangulation
    CGAL::draw(T);
  }

  return EXIT_SUCCESS;
}
