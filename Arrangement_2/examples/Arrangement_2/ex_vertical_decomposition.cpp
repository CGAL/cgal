//! \file examples/Arrangement_2/ex_vertical_decomposition.cpp
// Performing vertical decomposition of an arrangement.

#include <CGAL/MP_Float.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_vertical_decomposition.h>
#include <list>

typedef CGAL::MP_Float                                  Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                   Arrangement_2;
typedef Arrangement_2::Vertex_const_handle              Vertex_const_handle;
typedef Arrangement_2::Halfedge_const_handle            Halfedge_const_handle;

typedef std::list<Vertex_const_handle>                  Vertex_list;
typedef CGAL::Unique_hash_map<Vertex_const_handle,
                              std::pair<CGAL::Object,
                                        CGAL::Object> > Vertical_map;

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
  Vertical_map     v_map;
  Vertex_list      vs;

  CGAL::decompose (arr, v_map, std::back_inserter(vs));

  // Print the results.
  Vertex_list::const_iterator   vit;
  Vertex_const_handle           curr_v;
  Vertex_const_handle           vh;
  Halfedge_const_handle         hh;

  for (vit = vs.begin(); vit != vs.end(); ++vit)
  {
    curr_v = *vit;
    std::cout << "Vertex (" << curr_v->point() << ") : ";

    std::cout << " feature below: ";
    if (CGAL::assign (hh, v_map[curr_v].first))
    {
      if (! hh->is_fictitious())
        std::cout << '[' << hh->curve() << ']';
      else
        std::cout << "NONE";
    }
    else if (CGAL::assign (vh, v_map[curr_v].first))
      std::cout << '(' << vh->point() << ')';
    else
      std::cout << "EMPTY";

    std::cout << "   feature above: ";
    if (CGAL::assign (hh, v_map[curr_v].second))
    {
      if (! hh->is_fictitious())
        std::cout << '[' << hh->curve() << ']' << std::endl;
      else
        std::cout << "NONE" << std::endl;
    }
    else if (CGAL::assign (vh, v_map[curr_v].second))
      std::cout << '(' << vh->point() << ')' << std::endl;
    else
      std::cout << "EMPTY" << std::endl;
  }

  return (0);
}

