#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Snap_rounding_traits_2.h>
#include <CGAL/Snap_rounding_2.h>
#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_standard_toolbar.h>
#include <CGAL/IO/Qt_help_window.h>
#include <CGAL/IO/Qt_widget_layer.h>
#include <CGAL/IO/pixmaps/demoicon.xpm> 

typedef CGAL::Quotient<CGAL::MP_Float>         Number_type;
typedef CGAL::Cartesian<Number_type>           Rep;
typedef CGAL::Snap_rounding_traits_2<Rep>      Sr_traits;
typedef Rep::Segment_2                         Segment_2;
typedef Rep::Point_2                           Point_2;
typedef std::list<Segment_2>                   Segment_list_2;
typedef std::list<Point_2>                     Polyline_2;
typedef std::list<Polyline_2>                  Polyline_list_2;
typedef CGAL::Iso_rectangle_2<Rep>             Iso_rectangle_2;

typedef std::list<Segment_2>                   Segment_2_list;
typedef Segment_2_list::const_iterator         Segment_2_list_const_iterator;
typedef Segment_2_list::iterator               Segment_2_list_iterator;
typedef std::list<Point_2>                     Point_2_list;
typedef Point_2_list::const_iterator           Point_2_list_const_iterator;
typedef std::list<std::list<Point_2> >         Polyline_2_list;
typedef Polyline_2_list::const_iterator        Polyline_2_list_const_iterator;

#define MIN_X 0
#define MIN_Y 0
#define MAX_X 10
#define MAX_Y 10
#define PRECISION 0.5
