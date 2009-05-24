//! \file examples/Arrangement_on_surface_2/vertical_decomposition.cpp
// Performing vertical decomposition of an arrangement.

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_vertical_decomposition_2.h>
#include <list>

typedef CGAL::MP_Float                                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_linear_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Vertex_const_handle              Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle            Halfedge_const_handle;
typedef Arrangement_2::Face_const_handle                Face_const_handle;

typedef std::pair<Vertex_const_handle, std::pair<CGAL::Object, CGAL::Object> >
                                                        Vert_decomp_entry;
typedef std::list<Vert_decomp_entry>                    Vert_decomp_list;

int main ()
{
  // Construct the arrangement.
  Arrangement_2    arr;

  insert_non_intersecting_curve(arr, Segment_2(Point_2(1, 1), Point_2(3, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(1, 1), Point_2(2, 2)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 2), Point_2(3, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 2), Point_2(5, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(3, 2), Point_2(5, 0)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(2, 3), Point_2(3, 3)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(0, 3), Point_2(6, 4)));
  insert_non_intersecting_curve(arr, Segment_2(Point_2(4, 4), Point_2(4, 5)));

  // Perform vertical ray-shooting from every vertex and locate the feature
  // that lie below it and the feature that lies above it.
  Vert_decomp_list  vd_list;

  CGAL::decompose (arr, std::back_inserter(vd_list));

  // Print the results.
  Vert_decomp_list::const_iterator       vd_iter;
  std::pair<CGAL::Object, CGAL::Object>  curr;
  Vertex_const_handle                    vh;
  Halfedge_const_handle                  hh;
  Face_const_handle                      fh;

  for (vd_iter = vd_list.begin(); vd_iter != vd_list.end(); ++vd_iter) {
    curr = vd_iter->second;
    std::cout << "Vertex (" << vd_iter->first->point() << ") : ";

    std::cout << " feature below: ";
    if (CGAL::assign (vh, curr.first))
          std::cout << '(' << vh->point() << ')';
    else if (CGAL::assign (hh, curr.first))
      if (!hh->is_fictitious())
        std::cout << '[' << hh->curve() << ']';
      else
        std::cout << "NONE";
    else
      std::cout << "EMPTY";

    std::cout << "   feature above: ";
    if (CGAL::assign (vh, curr.second))
          std::cout << '(' << vh->point() << ')' << std::endl;
    else if (CGAL::assign (hh, curr.second))
      if (!hh->is_fictitious())
        std::cout << '[' << hh->curve() << ']' << std::endl;
      else
        std::cout << "NONE" << std::endl;
    else
      std::cout << "EMPTY" << std::endl;
  }

  return 0;
}
