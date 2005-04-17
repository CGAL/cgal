// An example program for using the Conic(CORE) traits.
// Should be used with CORE version 1.6x (the latest one).

#include <CGAL/Cartesian.h>
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Arr_conic_traits_2_core.h>
#include <CGAL/Arr_curve_data_traits_2.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Arr_2_default_dcel.h>
#include <CGAL/Arrangement_2.h>

typedef CORE::BigInt                                CfNT;
typedef CORE::Expr                                  CoNT;
typedef CGAL::Cartesian<CoNT>                       Kernel;
typedef CGAL::Arr_conic_traits_2<CfNT,Kernel>       Base_traits;
typedef Base_traits::Curve_2                        Base_curve_2;

typedef CGAL::Arr_curve_data_traits_2<Base_traits, int>   Traits;
typedef Traits::Point_2                                   Point_2;
typedef Traits::Curve_2                                   Curve_2;

typedef CGAL::Pm_default_dcel<Traits>               Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>             Pm_2;
typedef CGAL::Planar_map_with_intersections_2<Pm_2> Pmwx_2;

int main()
{
  // Create the arrangement using the "planar map with intersections" template.
  Pmwx_2                arr;  

  const int                  max_curves = 20;
  std::vector<Base_curve_2>  base_curves (max_curves);
  std::list<Curve_2>         curves;
  Base_curve_2               cv;
  int                        ind = 0;

  // Creation of a line segment with integer (of type CfNT) endpoints:
  // In this example we create the segment [(1, 1); (0, -3)].
  cv = Base_curve_2 (1, 1,             // Source coordinates.
		     0, -3);           // Target coordinates.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Creation of a line segment with "real" coordinates (of type CoNT).
  // In this case we supply the equation of the underlying line which is
  //   0*x^2 + 0*y^2 + 0*xy + a*x + b*y + c = 0,
  // where a, b, c must be integers (of type CfNT).
  // In this example we create the segment [(-0.5, -0.5); (sqrt(2), sqrt(2))],
  // which lies on the line y = x.
  CoNT                minus_half ("-0.5");
  CoNT                sqrt_2 = CGAL::sqrt(CoNT(2));
  Point_2             ps_1 (minus_half, minus_half);
  Point_2             pt_1 (sqrt_2, sqrt_2);

  cv = Base_curve_2 (0, 0, 0,          // The conic equation: in case of
		     1, -1, 0,         // a line, the first three coefficients
                                       // are always 0.
		     CGAL::COLLINEAR,  // The orientation.
		     ps_1,             // The source point.
		     pt_1);            // The target point.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Insert the circle centered at (1,1) with radius 1.
  // Notice that the center coordinates and the radius must always have
  // integer values (of type CfNT).
  cv = Base_curve_2 (1, 1,             // Center coordinates.
		     1);               // The radius.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Insert a circular arc that lies on a circle centered at (1,4)
  // with radius 3, going clockwise from (3, 4 - sqrt(5) to (-1, 4 - sqrt(5)).
  // Notice that the center coordinates and the radius must always have
  // integer values (of type CfNT), but the endpoints can have "real"
  // coordinates (of type CoNT).
  CoNT                y_coor = 4 - CGAL::sqrt(CoNT(5));
  Point_2             ps_2 (3, y_coor);
  Point_2             pt_2 (-1, y_coor);

  cv = Base_curve_2 (1, 4,             // Center coordinates.
		     3,                // The radius.
		     CGAL::CLOCKWISE,  // The orientation.
		     ps_2,             // The source point.
		     pt_2);            // The target point.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Insert a circle centered at (0.5,3) with radius sqrt(2).
  // The equation of this circle is:
  //   4*x^2 + 4*y^2 - 4*x - 24*y + 29 = 0
  // It is convenient to use this construction method when the circle center
  // has rational coordinates and its sqaured radius is rational.
  cv = Base_curve_2 (4, 4, 0, -4, -24, 29);
                                       // The 0 is the coefficient of xy
                                       // (always 0 for a circle).
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Insert an arc that lies on the circle centered at (-0.5,3) with radius 
  // sqrt(2). The equation of this circle is:
  //   4*x^2 + 4*y^2 + 4*x - 24*y + 29 = 0
  // The arc goes in a clockwise direction from (-0.5, 3 + sqrt(2)) to
  // (-0.5 - sqrt(2), 3) (from "12 o'clock" to "9 o'clock"). 
  Point_2             ps_3 (minus_half, 3 + sqrt_2);
  Point_2             pt_3 (minus_half - sqrt_2, 3);

  cv = Base_curve_2 (4, 4, 0, 4, -24, 29, // The curve coefficients.
		     CGAL::CLOCKWISE,     // The orientation.
		     ps_3,                // The source point.
		     pt_3);               // The target point.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;

  // Insert a circular arc, defined by the three points (0,-2), (2,-1)
  // and (1,1).
  // In this case, the three points must have integer coordinates (type CfNT).
  cv = Base_curve_2 (0, -2,               // Source coordinates.
		     2, -1,               // Mid-point coordinates.  
		     1, 1,                // Target coordinates.
		     true);               // Dummy parameter.
  base_curves[ind] = cv;
  curves.push_back (Curve_2 (cv, ind));
  ind++;
		    
  // Insert all curves to the arrangement.
  arr.insert (curves.begin(), curves.end());

  // Print out the number of vertices, edges and faces in the arrangement.
  std::cout << "Number of vertices: " 
	    << arr.number_of_vertices() << std::endl;
  
  std::cout << "Number of edges: " 
	    << arr.number_of_halfedges()/2 << std::endl;
  std::cout << "Number of faces: " 
	    << arr.number_of_faces() << std::endl;

  // Go over all vertices and for each vertex print the ID numbers of the
  // base curves that go through it.
  Pmwx_2::Vertex_iterator   vit;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); vit++)
  {
    // Scan the half-edges that has the current vertex as their endpoint.
    std::set<int>         indices;

    Pmwx_2::Halfedge_around_vertex_circulator 
      eit, first = (*vit).incident_halfedges();

    eit = first;
    do 
    {
      ind = (*eit).curve().get_data();

      // Keep track of IDs we haven't seen before.
      if (indices.find(ind) == indices.end())
	indices.insert (ind);

      eit++;

    } while (eit != first);

    // Disregard vertices with only one curve going through them - these are
    // not intersection points.
    if (indices.size() == 1)
      continue;

    // Print the vertex.
    const Point_2& p = (*vit).point();
    std::cout << "(" << CGAL::to_double(p.x()) 
	      << "," << CGAL::to_double(p.y()) << ") : ";

    // Print the indices.
    std::set<int>::const_iterator  ii;

    for (ii = indices.begin(); ii != indices.end(); ii++)
      std::cout << *ii << "  ";
    std::cout << std::endl;
  }

  return 0;
}
