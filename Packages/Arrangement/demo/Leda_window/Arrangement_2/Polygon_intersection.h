#include <CGAL/basic.h> 
#include <CGAL/leda_rational.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Polygon_2.h>

#include <CGAL/Arr_2_bases.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Pm_naive_point_location.h>

#include <list>

/////////////////////////////////////////////////////////////////////////
// Type definitions
/////////////////////////////////////////////////////////////////////////

typedef leda_rational                   NT;
typedef CGAL::Cartesian<NT>             R;

typedef CGAL::Arr_segment_traits_2<R>   Traits;
typedef Traits::Curve_2                 Curve;
typedef Traits::X_monotone_curve_2      X_curve;
typedef Traits::Point_2                 Point;
typedef CGAL::Arr_base_node<Curve, X_curve> Base_node;


typedef CGAL::Polygon_traits_2<R>                           Polygon_traits;
typedef std::list<Point>                                    Polygon_Container;
typedef CGAL::Polygon_2<Polygon_traits, Polygon_Container > Polygon;
typedef std::list<Polygon>                                  Polygon_list;

//a face with a counter for overlapping polygons (counter is initialized to -1)
struct Face_with_counter : public CGAL::Arr_2_face_base {
  Face_with_counter() : CGAL::Arr_2_face_base(), counter(-1) {}
  int counter;
};

//a DCEL with Face_with_counter
typedef CGAL::Pm_dcel<CGAL::Arr_2_vertex_base<Point>,
  CGAL::Arr_2_halfedge_base<Base_node >,
  Face_with_counter >  Dcel;

typedef CGAL::Arrangement_2<Dcel,Traits,Base_node >      Arr_2;
typedef CGAL::Pm_naive_point_location<Arr_2::Planar_map> Arr_naive_pl;

typedef Arr_2::Ccb_halfedge_circulator   Ccb_halfedge_circulator;
typedef Arr_2::Face_handle               Face_handle;
typedef Arr_2::Face_iterator             Face_iterator;
typedef Arr_2::Holes_iterator            Holes_iterator;
typedef std::list<Arr_2::Curve_iterator> ArrCurvesList;

/////////////////////////////////////////////////////////////////////////
// Function Prototypes
/////////////////////////////////////////////////////////////////////////

bool intersect_polygons(Polygon_list &in_poly_list,
                        Polygon_list &out_poly_list);
