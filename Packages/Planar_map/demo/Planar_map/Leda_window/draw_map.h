#ifndef __DRAW_MAP_H
#define __DRAW_MAP_H

#include <iostream>

// if LEDA is not installed, a message will be issued in runtime by demo.C.
#ifdef CGAL_USE_LEDA

#include "configuration"

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>  
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Pm_segment_traits_2.h>
#ifdef USE_RATIONAL
#include <CGAL/leda_rational.h>
#else
#if defined (USE_LEDA_RAT_KERNEL)
#include <CEP/Leda_rat_kernel/leda_rat_kernel_traits.h>
#endif
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/IO/Window_stream.h>

//#define CGAL_PM_DEBUG
#ifdef CGAL_PM_DEBUG
#include <CGAL/IO/Planar_map_iostream.h>
#endif

#ifdef CGAL_PM_TIMER
#include <CGAL/Timer.h>
#endif

#define BUNDLE 100
#define WIDE_PRECISION 10

#if defined(USE_RATIONAL) || defined(USE_LEDA_RAT_KERNEL)
#if defined(USE_RATIONAL) && defined(USE_LEDA_RAT_KERNEL)
#error only one kernel should be defined
#endif
typedef leda_rational                           number_type;
#else
typedef double                                  number_type; 
#endif

#ifdef USE_LEDA_RAT_KERNEL
typedef CGAL::leda_rat_kernel_traits            Rep;
#else
typedef CGAL::Cartesian<number_type>            Rep;
#endif

typedef CGAL::Pm_segment_traits_2<Rep>          Traits;

typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef Planar_map::Traits_wrap                 Traits_wrap;

typedef Planar_map::Vertex                      Vertex;
typedef Planar_map::Halfedge                    Halfedge;
typedef Planar_map::Face                        Face;

typedef Planar_map::Vertex_handle               Vertex_handle;
typedef Planar_map::Halfedge_handle             Halfedge_handle;
typedef Planar_map::Face_handle                 Face_handle;

typedef Planar_map::Vertex_iterator             Vertex_iterator;
typedef Planar_map::Halfedge_iterator           Halfedge_iterator;

typedef Planar_map::Ccb_halfedge_circulator     Ccb_halfedge_circulator;

typedef Traits::Point                           Pm_point;
typedef Traits::X_curve                         Pm_curve; 


extern  int draw_pm (Planar_map & pm , CGAL::Window_stream & W);

extern  bool Init (char *filename , Planar_map & pm) ;

extern  void win_border( double &x0 , double &x1 , double &y0 ,
                            Planar_map &pm);

extern  CGAL::Window_stream& operator<<(CGAL::Window_stream& os,
                                          Planar_map &M);

extern  void window_input(Planar_map & M, CGAL::Window_stream &W );

#ifdef CGAL_PM_TIMER
extern CGAL::Timer t_total,t_construction,t_insert,t_remove,t_locate,t_vertical;
extern int n_total,n_insert,n_remove,n_locate,n_vertical;
#endif

/* move to Eyals leda_rat ? */
#ifdef USE_LEDA_RAT_KERNEL
inline CGAL::Window_stream& operator<<(CGAL::Window_stream& os, const Pm_point& p){
    return os << leda_point(p.xcoordD(),p.ycoordD()); 
  } 
inline CGAL::Window_stream& operator<<(CGAL::Window_stream& os, const Pm_curve& c){
    leda_segment s(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D()); 
    return os << s; 
  }
#endif

#endif // CGAL_USE_LEDA

#endif // __DRAW_MAP_H

