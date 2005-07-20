#ifndef CGAL_TYPES_HEADER
#define CGAL_TYPES_HEADER

#include <CGAL/basic.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_help_window.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_curve_data_traits_2.h>
#include <CGAL/Arr_curve_origin_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Pm_with_intersections.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/IO/write_pm.h> 
#include <CGAL/IO/Pm_iostream.h> 
#include <CORE/BigInt.h>
#include <CGAL/CORE_Expr.h>
#include <CGAL/Planar_map_2.h> 
#include <CGAL/Bbox_2.h>
#include <CGAL/IO/Pm_drawer.h> 
#include <CGAL/IO/draw_pm.h> 

#include <CGAL/Pm_trapezoid_ric_point_location.h>
#include <CGAL/Pm_naive_point_location.h>
#include <CGAL/Pm_simple_point_location.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
//#include <CGAL/Pm_lenmarks_point_location.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/IO/Postscript_file_stream.h> 
#endif

#include <qcolordialog.h> 
#include <qcolor.h>  // color of faces (stored in curve data)
#include <iostream>

#include<CGAL/Gmpz.h>
#include<CGAL/Gmpq.h>

#include<vector>


enum TraitsType { SEGMENT_TRAITS, POLYLINE_TRAITS , CONIC_TRAITS};
enum SnapMode   { NONE , GRID , POINT};
enum Mode       { INSERT , DELETE , POINT_LOCATION , RAY_SHOOTING ,
                  DRAG , MERGE , SPLIT,FILLFACE};
enum ConicType  { CIRCLE , SEGMENT ,ELLIPSE , PARABOLA , HYPERBOLA};
enum Strategy   { NAIVE , SIMPLE , TRAP , WALK };

// default background color
const QColor def_bg_color(0,0,0);

// Coordinate related typedef - using inexact number type
typedef float                                              Coord_type;
typedef CGAL::Cartesian<Coord_type>                        Coord_kernel;
typedef Coord_kernel::Point_2                              Coord_point;
typedef Coord_kernel::Segment_2                            Coord_segment;
typedef Coord_kernel::Circle_2                             Coord_circle;


typedef CGAL::Polygon_2< Coord_kernel> Polygon;  // polygon is usefull for filling faces

// For the conic traits:
//#define COORD_SCALE  10

// Planar map typedef - using rational exact number type

//typedef CGAL::Quotient<CGAL::MP_Float>                     NT;
//typedef CGAL::Cartesian<NT>                                Kernel;


typedef CGAL::Gmpq                                         NT;
typedef CGAL::Cartesian<NT>                                Kernel;


typedef CORE::BigInt                                       CfNT;
typedef CGAL::Cartesian<CfNT>                              Int_kernel;
typedef CORE::Expr                                         CoNT;
typedef CGAL::Cartesian<CoNT>                              Alg_kernel;


template <class Info>
class Face_with_info : public CGAL::Pm_face_base {
  Info data;
  bool _visited;
  int _index;
  /*! array of info's for overlay issues */
  std::vector<Info> overlay_info;
public:

  // iterator for overlay_info vector
  typedef typename std::vector<Info>::iterator OverlayInfoIterator;

  Face_with_info() : CGAL::Pm_face_base(), data(),_visited(false),
                     _index(0) , overlay_info(20) {}

  Info info() { return data; }
  void set_info(Info i) { data = i; }
  bool visited() { return _visited; }
  void set_visited(bool b) { _visited = b; }
  int index() { return _index; }
  void set_index(int index) { _index = index; }
  Info get_overlay_info(int index) { return overlay_info[index]; }
  void set_overlay_info(int index, Info info) { overlay_info[index] = info; }
  void assign_overlay_info(std::vector<Info> info){ overlay_info = info; }

  OverlayInfoIterator OverlayInfoBegin() { return overlay_info.begin(); }
  OverlayInfoIterator OverlayInfoEnd() { return overlay_info.end(); }
};

template <class Traits , class Info>
class Dcel : public CGAL::Pm_dcel<
  CGAL::Pm_vertex_base<typename Traits::Point>,
  CGAL::Pm_halfedge_base<typename Traits::X_curve>,
  Face_with_info < Info >
> 
{
public:  // CREATION
  
  Dcel() {}
};

// forward decleration 
class Curve_data;

// Segments: 
typedef CGAL::Arr_segment_cached_traits_2<Kernel>       Base_seg_traits; 
typedef Base_seg_traits::Curve_2                        Pm_base_seg_2; 

typedef CGAL::Arr_curve_data_traits_2<Base_seg_traits,
                                      Curve_data >      Seg_traits; 
typedef Seg_traits::Curve_2                             Pm_seg_2; 
typedef Seg_traits::X_monotone_curve_2                  Pm_xseg_2;
typedef Seg_traits::Point_2                             Pm_seg_point_2;
typedef Dcel<Seg_traits , QColor>                       Seg_dcel;
typedef CGAL::Planar_map_2<Seg_dcel, Seg_traits>        Seg_pm;
typedef CGAL::Planar_map_with_intersections_2<Seg_pm>   Seg_arr;
typedef Seg_arr::Halfedge                               Seg_halfedge;
typedef Seg_arr::Locate_type                            Seg_locate_type;
typedef Seg_arr::Halfedge_handle                        Seg_halfedge_handle;
typedef Seg_arr::Face_handle                            Seg_face_handle;
typedef Seg_arr::Ccb_halfedge_circulator
  Seg_ccb_halfedge_circulator;
typedef Seg_arr::Holes_iterator                         Seg_holes_iterator;
typedef Seg_arr::Face_iterator                          Seg_face_iterator;
typedef std::list<Pm_seg_2*>                            Pm_seg_list;
typedef Pm_seg_list::const_iterator                     Pm_seg_const_iter;
typedef Pm_seg_list::iterator                           Pm_seg_iter;
typedef Seg_arr::Pmwx_change_notification
  Seg_arr_change_notification;

//point location
typedef CGAL::Pm_trapezoid_ric_point_location<Seg_pm>
  Seg_trap_point_location;
typedef CGAL::Pm_naive_point_location<Seg_pm>
  Seg_naive_point_location;
typedef CGAL::Pm_simple_point_location<Seg_pm>
  Seg_simple_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Seg_pm>
  Seg_walk_point_location;
//typedef CGAL::Pm_nearest_neighbor<Seg_pm>
// Seg_nearest_neighbor;
//typedef CGAL::Pm_lenmarks_point_location<Seg_pm,Nearest_neighbor>
//  Seg_lenmarks_point_location;

class Curve_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  Seg_halfedge_handle halfedge_handle;
  union Pointer { 
    Pm_base_seg_2 * m_curve; 
    Pm_xseg_2  *m_x_motonote_curve; 
  } m_ptr; 
}; 

// Polyline

//forward decleration
class Curve_pol_data;

typedef CGAL::Arr_polyline_traits_2<Base_seg_traits>    Base_pol_traits;
typedef Base_pol_traits::Curve_2                        Pm_base_pol_2;

typedef CGAL::Arr_curve_data_traits_2<Base_pol_traits,
                                      Curve_pol_data>   Pol_traits;
typedef Pol_traits::Curve_2                             Pm_pol_2;
typedef Pol_traits::X_monotone_curve_2                  Pm_xpol_2;

typedef Pol_traits::Point_2                             Pm_pol_point_2;
typedef Dcel<Pol_traits,QColor>                         Pol_dcel;
typedef CGAL::Planar_map_2<Pol_dcel, Pol_traits>        Pol_pm;
typedef CGAL::Planar_map_with_intersections_2<Pol_pm>   Pol_arr;
typedef Pol_arr::Locate_type                            Pol_locate_type;
typedef Pol_arr::Halfedge_handle                        Pol_halfedge_handle;
typedef Pol_arr::Face_handle                            Pol_face_handle;
typedef Pol_arr::Ccb_halfedge_circulator
  Pol_ccb_halfedge_circulator;
typedef Pol_arr::Holes_iterator                         Pol_holes_iterator;
typedef Pol_arr::Halfedge                               Pol_halfedge;
typedef Pol_arr::Face_iterator                          Pol_face_iterator;
typedef Pol_arr::Pmwx_change_notification
  Pol_arr_change_notification;

typedef std::list<Pm_pol_2*>                            Pm_pol_list;
typedef Pm_pol_list::const_iterator                     Pm_pol_const_iter;
typedef Pm_pol_list::iterator                           Pm_pol_iter;

//point location
typedef CGAL::Pm_trapezoid_ric_point_location<Pol_pm>
  Pol_trap_point_location;
typedef CGAL::Pm_naive_point_location<Pol_pm>
  Pol_naive_point_location;
typedef CGAL::Pm_simple_point_location<Pol_pm>
  Pol_simple_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Pol_pm>
  Pol_walk_point_location;
//typedef CGAL::Pm_nearest_neighbor<Pol_pm>
// Pol_nearest_neighbor;
//typedef CGAL::Pm_lenmarks_point_location<Pol_pm,Nearest_neighbor>
// Pol_lenmarks_point_location;

class Curve_pol_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  Pol_halfedge_handle halfedge_handle;
  union Pointer { 
    Pm_base_pol_2 * m_curve; 
    Pm_xpol_2* m_x_motonote_curve; 
  } m_ptr;
};

// Conics

//forward decleration
class Curve_conic_data;


typedef CGAL::Arr_conic_traits_2<Int_kernel, Alg_kernel>   Base_conic_traits;
typedef Base_conic_traits::Curve_2                         Pm_base_conic_2;
typedef Base_conic_traits::Int_point_2                     Int_point_2;
typedef Base_conic_traits::Int_segment_2                   Int_segment_2;
typedef Base_conic_traits::Int_circle_2                    Int_circle_2;
typedef Base_conic_traits::Int_line_2                      Int_line_2;

typedef CGAL::Arr_curve_data_traits_2<Base_conic_traits,
                                      Curve_conic_data> Conic_traits;
typedef Conic_traits::Curve_2                           Pm_conic_2;
typedef Conic_traits::X_monotone_curve_2                Pm_xconic_2;
typedef Conic_traits::Point_2                           Pm_conic_point_2;
typedef Dcel<Conic_traits,QColor>                       Conic_dcel;
typedef CGAL::Planar_map_2<Conic_dcel, Conic_traits>    Conic_pm;
typedef CGAL::Planar_map_with_intersections_2<Conic_pm> Conic_arr;
typedef Conic_arr::Locate_type                          Conic_locate_type;
typedef Conic_arr::Halfedge_handle                      Conic_halfedge_handle;
typedef Conic_arr::Face_handle                          Conic_face_handle;
typedef Conic_arr::Ccb_halfedge_circulator
  Conic_ccb_halfedge_circulator;
typedef Conic_arr::Holes_iterator                       Conic_holes_iterator;
typedef CGAL::Pm_file_scanner<Conic_arr>                Pm_scanner; 
typedef Conic_arr::Halfedge                             Conic_halfedge;
typedef Conic_arr::Face_iterator                        Conic_face_iterator;
typedef Conic_arr::Pmwx_change_notification
  Conic_arr_change_notification;

typedef std::list<Pm_xconic_2*>                         Pm_xconic_list;
typedef Pm_xconic_list::const_iterator                  Pm_xconic_const_iter;
typedef Pm_xconic_list::iterator                        Pm_xconic_iter;

//point location
typedef CGAL::Pm_trapezoid_ric_point_location<Conic_pm>
  Conic_trap_point_location;
typedef CGAL::Pm_naive_point_location<Conic_pm>
  Conic_naive_point_location;
typedef CGAL::Pm_simple_point_location<Conic_pm>
  Conic_simple_point_location;
typedef CGAL::Pm_walk_along_line_point_location<Conic_pm>
  Conic_walk_point_location;
//typedef CGAL::Pm_nearest_neighbor<Conic_pm>
// Conic_nearest_neighbor;
//typedef CGAL::Pm_lenmarks_point_location<Conic_pm,Nearest_neighbor>
// Conic_lenmarks_point_location;


class Curve_conic_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  ConicType m_ct; 
  Conic_halfedge_handle halfedge_handle;
  union Pointer { 
    Pm_base_conic_2 * m_curve; 
    Pm_xconic_2 * m_x_motonote_curve; 
  } m_ptr; 
};

#endif
