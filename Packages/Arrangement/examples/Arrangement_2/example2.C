// examples/Arrangement_2/example2.C
// ---------------------------------

#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/circulator.h>

typedef CGAL::Quotient<CGAL::MP_Float>                  NT;
typedef CGAL::Cartesian<NT>                             Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits;
typedef Traits::Point_2                                 Point_2;
typedef Traits::Curve_2                                 Curve_2;
typedef Traits::X_monotone_curve_2                      X_monotone_curve_2;
typedef CGAL::Arr_2_default_dcel<Traits>                Dcel;
typedef CGAL::Arrangement_2<Dcel,Traits>                Arr_2;

int main() 
{
  Arr_2 arr;

  // Insertion of the curves
  arr.insert(Curve_2(Point_2(0, 0), Point_2(2, 2)));
  arr.insert(Curve_2(Point_2(1, 1), Point_2(3, 3)));

  // Traversal of the halfedges
  Arr_2::Halfedge_const_iterator hit = arr.halfedges_begin(),
      hit_end = arr.halfedges_end();

  for (; hit != hit_end; ++hit, ++hit)
  {
    // We skip the adjacent twin halfedge
    Arr_2::Overlap_const_circulator occ = hit->overlap_edges(), occ_end = occ;
    int count = 0;
    CGAL_For_all(occ, occ_end) { ++count; }

    if (count == 1) 
      std::cout << "Edge " << occ->x_curve() << " is covered by a single edge."
                << std::endl;
    else
      std::cout << "Edge " << occ->x_curve() << " is covered by " << count
                << " edges." << std::endl;
  }

  return 0;
}
