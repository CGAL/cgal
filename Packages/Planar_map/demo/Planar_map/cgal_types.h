#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/squared_distance_2.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Pm_segment_traits_2.h>
#include <CGAL/Planar_map_2.h>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Pm_segment_traits_2<Kernel>       Traits;

typedef Traits::Point_2                         Point_2;
typedef Traits::X_monotone_curve_2              Curve;

typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;

typedef Planar_map::Locate_type                 Locate_type;

typedef Planar_map::Halfedge_handle             Halfedge_handle;
typedef Planar_map::Face_handle                 Face_handle;

typedef Planar_map::Halfedge_iterator           Halfedge_iterator;
typedef Planar_map::Holes_iterator              Holes_iterator;
typedef Planar_map::Ccb_halfedge_circulator     Ccb_halfedge_circulator;
