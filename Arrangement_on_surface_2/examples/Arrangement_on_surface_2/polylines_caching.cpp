//! \file examples/Arrangement_on_surface_2/polylines_caching.cpp
// Constructing an arrangement of polylines with caching traits

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_caching_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Polygon_2.h>
#include <boost/function_output_iterator.hpp>
#include <vector>
#include <list>

#include "arr_print.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_2 = Kernel::Point_2;
using Line_2 = Kernel::Line_2;
using Point_2_vector = std::vector<Point_2>;
using Polygon_2 = CGAL::Polygon_2<Kernel>;

int main()
{
  // Using a vector of points
  {
    using Geom_traits_2 = CGAL::Arr_caching_polyline_traits_2<Kernel, Point_2_vector>;
    using Arrangement_2 = CGAL::Arrangement_2<Geom_traits_2>;
    using Curve_2 = typename Geom_traits_2::Curve_2;
    using X_monotone_curve_2 = typename Geom_traits_2::X_monotone_curve_2;

    Geom_traits_2 traits;
    Arrangement_2 arr(&traits);

    Point_2_vector points;
    points.push_back(Point_2(1, 3));
    points.push_back(Point_2(0, 2));
    points.push_back(Point_2(1, 0));
    points.push_back(Point_2(2, 1));
    points.push_back(Point_2(3, 0));
    points.push_back(Point_2(4, 1));
    points.push_back(Point_2(5, 0));
    points.push_back(Point_2(6, 2));
    points.push_back(Point_2(5, 3));
    points.push_back(Point_2(4, 2));

    auto construct_curve_2 = traits.construct_curve_2_object();
    Curve_2 curve = construct_curve_2(points);

    auto make_x_monotone_2 = traits.make_x_monotone_2_object();
    std::vector<X_monotone_curve_2> x_monotone_curves;
    make_x_monotone_2
      (curve,
       boost::make_function_output_iterator
       ([&](const CGAL::Object& obj)
        {
          const X_monotone_curve_2* xc = CGAL::object_cast<X_monotone_curve_2>(&obj);
          CGAL_assertion(xc);
          x_monotone_curves.push_back(*xc);
        }));

    for (const X_monotone_curve_2& xc : x_monotone_curves)
      insert (arr, xc);
  }

  // Using a Polygon_2
  {
    using Geom_traits_2 = CGAL::Arr_caching_polyline_traits_2<Kernel, Polygon_2>;
    using Arrangement_2 = CGAL::Arrangement_2<Geom_traits_2>;

    using Curve_2 = typename Geom_traits_2::Curve_2;
    using X_monotone_curve_2 = typename Geom_traits_2::X_monotone_curve_2;

    Geom_traits_2 traits;
    Arrangement_2 arr(&traits);

    Polygon_2 polygon;
    polygon.push_back (Point_2(0, 0));
    polygon.push_back (Point_2(2, 4));
    polygon.push_back (Point_2(3, 0));
    polygon.push_back (Point_2(4, 4));
    polygon.push_back (Point_2(6, 0));

    auto construct_curve_2 = traits.construct_curve_2_object();
    Curve_2 curve = construct_curve_2(polygon);

    auto make_x_monotone_2 = traits.make_x_monotone_2_object();
    std::vector<X_monotone_curve_2> x_monotone_curves;
    make_x_monotone_2
      (curve,
       boost::make_function_output_iterator
       ([&](const CGAL::Object& obj)
        {
          const X_monotone_curve_2* xc = CGAL::object_cast<X_monotone_curve_2>(&obj);
          CGAL_assertion(xc);
          x_monotone_curves.push_back(*xc);
        }));

    for (const X_monotone_curve_2& xc : x_monotone_curves)
      insert (arr, xc);
  }

  return EXIT_SUCCESS;
}
