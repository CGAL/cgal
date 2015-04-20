//! \file examples/Arrangement_on_surface_2/polylines.cpp
// Constructing an arrangement of polylines.

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <vector>
#include <list>

/*
  Define the Arrangement traits class to be used. You can either use some user
  defined kernel and Segment_traits_2 or the defaults.
 */

// Instantiate the traits class using a user-defined kernel
// and Segment_traits_2.
typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Segment_traits_2;
typedef CGAL::Arr_polyline_traits_2<Segment_traits_2>     Geom_traits_2;

// Identical instantiation can be achieved using the default Kernel:
// typedef CGAL::Arr_polyline_traits_2<>                    Geom_traits_2;


typedef Geom_traits_2::Point_2                            Point_2;
typedef Geom_traits_2::Segment_2                          Segment_2;
typedef Geom_traits_2::Curve_2                            Polyline_2;
typedef CGAL::Arrangement_2<Geom_traits_2>                Arrangement_2;
typedef Geom_traits_2::X_monotone_curve_2                 X_monotone_polyline;
typedef Geom_traits_2::X_monotone_subcurve_2              X_monotone_subcurve;

int main(int argc, char* argv[])
{
  Geom_traits_2 traits;
  Arrangement_2 arr(&traits);

  Geom_traits_2::Construct_x_monotone_curve_2 x_mono_polyline_construct =
    traits.construct_x_monotone_curve_2_object();

  std::vector<X_monotone_subcurve> x_segments;
  x_segments.push_back(X_monotone_subcurve( Point_2(0,0), Point_2(1,1) ));
  x_segments.push_back(X_monotone_subcurve( Point_2(1,1), Point_2(10,10) ));
  x_segments.push_back(X_monotone_subcurve( Point_2(10,10), Point_2(15,20) ));
  X_monotone_polyline pi1 = x_mono_polyline_construct(x_segments.begin(),
                                                      x_segments.end() );

  std::cout << "polline is: " << pi1 << std::endl;

  Point_2 src(atoi(argv[1]), atoi(argv[2]));
  Point_2 tgt(atoi(argv[3]), atoi(argv[4]));
  X_monotone_polyline trimmed_polyline = traits.trim_2_object()(pi1, src, tgt);
  std::cout << "trimmed polline is: " << trimmed_polyline << std::endl;

  return 0;
}
