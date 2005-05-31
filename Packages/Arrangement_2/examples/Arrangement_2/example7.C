// file: examples/Arrangement_2/example1.C


//#include "short_names.h"

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_observer.h>

enum Segment_color
{
  RED,
  BLUE
};

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;
typedef CGAL::Cartesian<Number_type>                  Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Segment_traits_2;
typedef Segment_traits_2::Curve_2                     Segment_2;
typedef CGAL::Arr_consolidated_curve_data_traits_2
                   <Segment_traits_2, Segment_color>  Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::Curve_2                             Colored_segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef CGAL::Arr_naive_point_location<Arrangement_2> Naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location
                                      <Arrangement_2> Walk_point_location;

int main ()
{
  // Construct the arrangement using the specialized insertion functions.
  Arrangement_2         arr;
  Naive_point_location  naive_pl (arr);

  Segment_2         seg1 (Point_2(-1, -1), Point_2(1, 3));
  Colored_segment_2 cs1 (seg1, RED);
  Segment_2         seg2 (Point_2(2, 0), Point_2(3, 3));
  Colored_segment_2 cs2 (seg2, RED);
  Segment_2         seg3 (Point_2(-1, 3), Point_2(4, 1));
  Colored_segment_2 cs3 (seg3, BLUE);
  Segment_2         seg4 (Point_2(-1, 0), Point_2(4, 1));
  Colored_segment_2 cs4 (seg4, BLUE);

  arr_insert (arr, naive_pl, cs1);
  arr_insert (arr, naive_pl, cs2);
  arr_insert (arr, naive_pl, cs3);
  arr_insert (arr, naive_pl, cs4);

  std::cout << "V = " << arr.number_of_vertices()
	    << ",  E = " << arr.number_of_edges() 
	    << ",  F = " << arr.number_of_faces() << std::endl;

  // Go over all vertices and print just the ones corresponding to intersection
  // points between RED segments and BLUE segments.
  Arrangement_2::Vertex_const_iterator   vit;
  Segment_color                          color;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    // Go over the incident edges of the current vertex and examine their
    // colors.
    bool       has_red = false;
    bool       has_blue = false;

    Arrangement_2::Halfedge_around_vertex_const_circulator  eit, first;

    eit = first = (*vit).incident_halfedges();
    do
    {
      // Get the color of the current half-edge.
      color = (*eit).curve().get_data();

      if (color == RED)
	has_red = true;
      else if (color == BLUE)
	has_blue = true;

      ++eit;
    } while (eit != first);

    // Print the vertex only if RED and BLUE half-edges were found.
    if (has_red && has_blue)
    {
      std::cout << "Red-blue intersection at (" << (*vit).point() << ")" 
		<< std::endl;
    }
  }

  // Copy the arrangement and color all edges blue.
  Arrangement_2                  arr2 = arr;
  Arrangement_2::Edge_iterator   eit;

  for (eit = arr2.edges_begin(); eit != arr2.edges_end(); ++eit)
  {
    if ((*eit).curve().get_data() == RED)
    {
      arr2.modify_edge (*eit,
			Traits_2::X_monotone_curve_2 ((*eit).curve(), BLUE));
    }
  }

  // Perform some ray-shooting queries.
  Point_2                              q = Point_2 (3, 2);
  Walk_point_location                  walk_pl (arr2);
  CGAL::Object                         obj;
  Arrangement_2::Vertex_const_handle   vh;
  Arrangement_2::Halfedge_const_handle hh;

  obj = walk_pl.ray_shoot_up (q);

  std::cout << "Shoot up from (" << q << ") : ";
  if (CGAL::assign (hh, obj))
  {
    std::cout << "Edge " << hh.curve() << std::endl;
  }
  else if (CGAL::assign (vh, obj))
  {
    std::cout << "Vertex " << vh.point() << std::endl;
  }
  else
  {
    std::cout << "Nothing." << std::endl;    
  }
   
  obj = walk_pl.ray_shoot_down (q);

  std::cout << "Shoot down from (" << q << ") : ";
  if (CGAL::assign (hh, obj))
  {
    std::cout << "Edge " << hh.curve() << std::endl;
  }
  else if (CGAL::assign (vh, obj))
  {
    std::cout << "Vertex " << vh.point() << std::endl;
  }
  else
  {
    std::cout << "Nothing." << std::endl;    
  }
 
  return (0);
}

