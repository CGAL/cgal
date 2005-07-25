// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

#ifndef CGAL_TYPES_HEADER
#define CGAL_TYPES_HEADER

#include <CGAL/basic.h>

#include <CGAL/IO/Qt_widget.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_help_window.h>

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_default_dcel.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_consolidated_curve_data_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Bbox_2.h>

#include <CGAL/Arr_trapezoid_ric_point_location.h>
#include <CGAL/Arr_naive_point_location.h>
#include <CGAL/Arr_walk_along_line_point_location.h>
#include <CGAL/Arr_landmarks_point_location.h>

#ifdef CGAL_USE_LEDA
#include <CGAL/IO/Postscript_file_stream.h> 
#endif

#include <qcolordialog.h> 
#include <qcolor.h>  // color of faces (stored in curve data)
#include <iostream>

#include <CGAL/Gmpz.h>
#include <CGAL/Gmpq.h>

#include <vector>

#include <CGAL/Arr_observer.h>


enum TraitsType { SEGMENT_TRAITS, POLYLINE_TRAITS , CONIC_TRAITS};
enum SnapMode   { NONE , GRID , POINT};
enum Mode       { INSERT , DELETE , POINT_LOCATION , RAY_SHOOTING_UP ,
                  RAY_SHOOTING_DOWN, DRAG , MERGE , SPLIT,FILLFACE};
enum ConicType  { CIRCLE , SEGMENT ,ELLIPSE , PARABOLA , HYPERBOLA};
enum Strategy   { NAIVE , TRAP , WALK, LANDMARKS };

// default background color
const QColor def_bg_color(0,0,0);

// Coordinate related typedef - using inexact number type
typedef float                                              Coord_type;
typedef CGAL::Cartesian<Coord_type>                        Coord_kernel;
typedef Coord_kernel::Point_2                              Coord_point;
typedef Coord_kernel::Segment_2                            Coord_segment;
typedef Coord_kernel::Circle_2                             Coord_circle;


typedef CGAL::Polygon_2<Coord_kernel> Polygon;  // polygon is usefull for filling faces

typedef CGAL::Gmpq                                         NT;
typedef CGAL::Cartesian<NT>                                Kernel;


template <class Info>
class Face_with_info : public CGAL::Arr_face_base {
  Info data;
  bool _visited;
  int _index;
  /*! array of info's for overlay issues */
  std::vector<Info> overlay_info;
public:

  // iterator for overlay_info vector
  typedef typename std::vector<Info>::iterator OverlayInfoIterator;

  Face_with_info() : CGAL::Arr_face_base(), data(),_visited(false),
                     _index(0) , overlay_info(20) {}

  Info info() const { return data; }
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
class Dcel : public CGAL::Arr_dcel<
  CGAL::Arr_vertex_base<typename Traits::Point_2>,
  CGAL::Arr_halfedge_base<typename Traits::X_monotone_curve_2>,
  Face_with_info < Info >
> 
{
public:  // CREATION
  
  Dcel() {}
};

// forward decleration 
class Curve_data;

// Segments: 
typedef CGAL::Arr_segment_traits_2<Kernel>              Base_seg_traits; 
typedef Base_seg_traits::Curve_2                        Arr_base_seg_2; 

typedef CGAL::Arr_consolidated_curve_data_traits_2<Base_seg_traits,
                                                   Curve_data > 
                                                        Seg_traits; 
typedef Seg_traits::Curve_2                             Arr_seg_2; 
typedef Seg_traits::X_monotone_curve_2                  Arr_xseg_2;
typedef Seg_traits::Point_2                             Arr_seg_point_2;
typedef Dcel<Seg_traits , QColor>                       Seg_dcel;
typedef CGAL::Arrangement_2<Seg_traits, Seg_dcel>       Seg_arr;
typedef Seg_arr::Halfedge                               Seg_halfedge;
typedef Seg_arr::Halfedge_handle                        Seg_halfedge_handle;
typedef Seg_arr::Face_handle                            Seg_face_handle;
typedef Seg_arr::Ccb_halfedge_circulator                
  Seg_ccb_halfedge_circulator;
typedef Seg_arr::Holes_iterator                         Seg_holes_iterator;
typedef Seg_arr::Face_iterator                          Seg_face_iterator;
typedef std::list<Arr_seg_2*>                           Arr_seg_list;
typedef Arr_seg_list::const_iterator                    Arr_seg_const_iter;
typedef Arr_seg_list::iterator                          Arr_seg_iter;


//point location
typedef CGAL::Arr_trapezoid_ric_point_location<Seg_arr>
  Seg_trap_point_location;
typedef CGAL::Arr_naive_point_location<Seg_arr>
  Seg_naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Seg_arr>
  Seg_walk_point_location;
typedef CGAL::Arr_landmarks_point_location<Seg_arr>
  Seg_lanmarks_point_location;



class Curve_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  union Pointer
  { 
    Arr_base_seg_2    *m_curve; 
    Arr_xseg_2        *m_x_motonote_curve; 
  }m_ptr; 

   bool operator== (const Curve_data& ccd) const
  {
    return (m_index == ccd.m_index && m_type == ccd.m_type);
  }
}; 

// Polyline

//forward decleration
class Curve_pol_data;

typedef CGAL::Arr_polyline_traits_2<Base_seg_traits>    Base_pol_traits;
typedef Base_pol_traits::Curve_2                        Arr_base_pol_2;

typedef CGAL::Arr_consolidated_curve_data_traits_2<Base_pol_traits,
                                                   Curve_pol_data>   
                                                        Pol_traits;
typedef Pol_traits::Curve_2                             Arr_pol_2;
typedef Pol_traits::X_monotone_curve_2                  Arr_xpol_2;

typedef Pol_traits::Point_2                             Arr_pol_point_2;
typedef Dcel<Pol_traits,QColor>                         Pol_dcel;
typedef CGAL::Arrangement_2<Pol_traits, Pol_dcel>       Pol_arr;
typedef Pol_arr::Halfedge_handle                        Pol_halfedge_handle;
typedef Pol_arr::Face_handle                            Pol_face_handle;
typedef Pol_arr::Ccb_halfedge_circulator
  Pol_ccb_halfedge_circulator;
typedef Pol_arr::Holes_iterator                         Pol_holes_iterator;
typedef Pol_arr::Halfedge                               Pol_halfedge;
typedef Pol_arr::Face_iterator                          Pol_face_iterator;

typedef std::list<Arr_pol_2*>                            Arr_pol_list;
typedef Arr_pol_list::const_iterator                     Arr_pol_const_iter;
typedef Arr_pol_list::iterator                           Arr_pol_iter;

//point location
typedef CGAL::Arr_trapezoid_ric_point_location<Pol_arr>
  Pol_trap_point_location;
typedef CGAL::Arr_naive_point_location<Pol_arr>
  Pol_naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Pol_arr>
  Pol_walk_point_location;
typedef CGAL::Arr_landmarks_point_location<Pol_arr>
  Pol_lanmarks_point_location;


class Curve_pol_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  Pol_halfedge_handle halfedge_handle;
  union Pointer { 
    Arr_base_pol_2 * m_curve; 
    Arr_xpol_2* m_x_motonote_curve; 
  } m_ptr;

   bool operator== (const Curve_pol_data& ccd) const
  {
    return (m_index == ccd.m_index && m_type == ccd.m_type &&
            halfedge_handle == ccd.halfedge_handle);
  }
};

// Conics

//forward decleration
class Curve_conic_data;

typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef CGAL::Arr_conic_traits_2<Rat_kernel, 
                                 Alg_kernel,
                                 Nt_traits>           Base_conic_traits;

typedef  Base_conic_traits::Curve_2                    Arr_base_conic_2;
typedef  Base_conic_traits::Rat_point_2                Rat_point_2;
typedef  Base_conic_traits::Rat_segment_2              Rat_segment_2;
typedef  Base_conic_traits::Rat_circle_2               Rat_circle_2;
typedef  Base_conic_traits::Rat_line_2                 Rat_line_2;

typedef CGAL::Arr_consolidated_curve_data_traits_2<Base_conic_traits,
                                      Curve_conic_data> Conic_traits;
typedef Conic_traits::Curve_2                           Arr_conic_2;
typedef Conic_traits::X_monotone_curve_2                Arr_xconic_2;
typedef Conic_traits::Point_2                           Arr_conic_point_2;
typedef Dcel<Conic_traits,QColor>                       Conic_dcel;
typedef CGAL::Arrangement_2<Conic_traits, Conic_dcel>   Conic_arr;
typedef Conic_arr::Halfedge_handle                      Conic_halfedge_handle;
typedef Conic_arr::Face_handle                          Conic_face_handle;
typedef Conic_arr::Ccb_halfedge_circulator
  Conic_ccb_halfedge_circulator;
typedef Conic_arr::Holes_iterator                       Conic_holes_iterator;
//typedef CGAL::Arr_file_scanner<Conic_arr>                Arr_scanner; 
typedef Conic_arr::Halfedge                             Conic_halfedge;
typedef Conic_arr::Face_iterator                        Conic_face_iterator;

typedef std::list<Arr_xconic_2*>                         Arr_xconic_list;
typedef Arr_xconic_list::const_iterator                  Arr_xconic_const_iter;
typedef Arr_xconic_list::iterator                        Arr_xconic_iter;

//point location
typedef CGAL::Arr_trapezoid_ric_point_location<Conic_arr>
  Conic_trap_point_location;
typedef CGAL::Arr_naive_point_location<Conic_arr>
  Conic_naive_point_location;
typedef CGAL::Arr_walk_along_line_point_location<Conic_arr>
  Conic_walk_point_location;
typedef CGAL::Arr_landmarks_point_location<Conic_arr>
 Conic_lanmarks_point_location;


class Curve_conic_data { 
public:
  enum Type {LEAF, INTERNAL}; 
  int m_index; 
  Type m_type; 
  ConicType m_ct; 
  Conic_halfedge_handle halfedge_handle;
  union Pointer { 
    Arr_base_conic_2 * m_curve; 
    Arr_xconic_2 * m_x_motonote_curve; 
  } m_ptr; 

  bool operator== (const Curve_conic_data& ccd) const
  {
    return (m_index == ccd.m_index && m_type == ccd.m_type &&
            m_ct == ccd.m_ct && halfedge_handle == ccd.halfedge_handle);
  }
};

template <class Arrangement_>
class My_observer : public CGAL::Arr_observer<Arrangement_>
{
public:

  typedef Arrangement_                        Arrangement;
  typedef CGAL::Arr_observer<Arrangement>     Arr_observer;
  typedef typename Arrangement::Face_handle   Face_handle;

  My_observer (Arrangement& arr) : Arr_observer (arr)
  {}

   virtual void after_split_face (Face_handle  f ,
                                  Face_handle  new_f ,
                                  bool         is_hole)
  {
    new_f ->set_info(f->info());
  }

};



#endif
