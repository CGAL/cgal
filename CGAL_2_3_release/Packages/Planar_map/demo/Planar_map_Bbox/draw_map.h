#ifndef __DRAW_MAP_H
#define __DRAW_MAP_H

// if LEDA is not installed, a message will be issued in runtime by demo.C.
#ifdef CGAL_USE_LEDA

#include "configuration"

#define BUNDLE 100
#define WIDE_PRECISION 10

#include <CGAL/Cartesian.h>
#include <CGAL/squared_distance_2.h>  
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Segment_2.h>

#ifdef USE_RATIONAL
#include <CGAL/leda_rational.h>
#endif
#include <CGAL/Pm_straight_exact_traits.h>
#include <CGAL/IO/Straight_2_stream.h>
#include <CGAL/IO/Pm_straight_exact_traits_stream.h>

#ifdef CGAL_PM_DEBUG
#include <CGAL/Pm_segment_exact_traits.h>
#include <CGAL/Pm_traits_checker.h>
#endif

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#ifdef CGAL_PM_DEBUG

#include <CGAL/IO/Planar_map_iostream.h>

#endif //CGAL_PM_DEBUG

#include <CGAL/IO/Window_stream.h>

#ifdef CGAL_PM_TIMER
#include <CGAL/Timer.h>
#endif

#ifdef USE_RATIONAL
typedef leda_rational                           number_type;
#else
typedef double                                  number_type;
#endif

typedef CGAL::Cartesian<number_type>            Rep;
typedef Rep::FT                                 FT;
typedef Rep::RT					RT;

typedef CGAL::Pm_straight_exact_traits<Rep>      Traits;
typedef CGAL::Pm_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>         Planar_map;
typedef Planar_map::Traits_wrap	Traits_wrap;

typedef CGAL::Window_stream			Window;

typedef Planar_map::Vertex                     Vertex;
typedef Planar_map::Halfedge                   Halfedge;
typedef Planar_map::Face                       Face;
typedef Planar_map::Vertex_handle              Vertex_handle;
typedef Planar_map::Halfedge_handle            Halfedge_handle;
typedef Planar_map::Face_handle                Face_handle;
typedef Planar_map::Locate_type		Locate_type;
typedef Planar_map::Vertex_iterator            Vertex_iterator;
typedef Planar_map::Halfedge_iterator          Halfedge_iterator;
typedef Planar_map::Ccb_halfedge_circulator    Ccb_halfedge_circulator;
typedef Planar_map::Point_location_base Point_location_base;
typedef Planar_map::Bounding_box_base Bounding_box_base;

typedef Traits::Point                          Point;
typedef Traits::X_curve                        X_curve; 
typedef Traits::X_unbounded_curve              Line;
typedef Traits::X_bounded_curve                Segment;
typedef Traits::X_target_unbounded_curve       Ray;

typedef Traits::X_curve_container X_curve_container;
typedef Traits::Point_container Point_container;
typedef X_curve_container::iterator X_curve_iterator;
typedef Point_container::iterator Point_iterator;

Window& read_line_ray_segment_plus(Window& w, const Point& s, Point& t, int& type,int& key);

/* Window& read_line(const Point& p , Line& l , Planar_map &m, Window& w);
Window& read_ray(const Point& p , Ray& s, Planar_map &pm, Window& w);
Window& read_segment(const Point& p , Segment& s, Planar_map &pm, Window& w); */
void find_face (Planar_map &pm,const Point& p , X_curve_container& l);
bool cooriented(const Traits_wrap * ,const Halfedge_handle& h,const X_curve& cv);
void split_edge(Planar_map &m, const Halfedge_handle& h,
				X_curve& cv1,X_curve& cv2,
				const Point& p);
Locate_type vertical_ray_shoot (Planar_map &m,const Point& p,Point & q, 
								bool up,Halfedge_handle& e);
void draw_arrow (CGAL::Window_stream & W ,const Point& p1, const Point& p2, bool black);
Vertex_handle find_closest_vertex(Planar_map &m, const Point& p);
void redraw(leda_window* wp, double x0, double y0, double x1, double y1); 

extern bool Init (char *filename , Planar_map & pm) ;
extern void win_border( double &x0 , double &x1 , double &y0 ,
						   Planar_map &pm);
//extern void window_input(Planar_map & M, CGAL::Window_stream &W );
extern int curve_type,action_type,display_mode,mbutton;
extern Planar_map* mp;
//,available;
extern	Point lastp,lastq;

#ifdef CGAL_PM_TIMER

extern CGAL::Timer t_total,t_construction,t_insert,t_remove,t_locate,t_vertical;
extern int n_total,n_insert,n_remove,n_locate,n_vertical;

#endif //CGAL_PM_TIMER

/* move to Eyal's leda_rat ?
#ifdef USE_LEDA_RAT_KERNEL

 inline CGAL::Window_stream& operator<<(CGAL::Window_stream& os, const Point& p){
 return os << leda_point(p.xcoordD(),p.ycoordD()); 
 } 
 inline CGAL::Window_stream& operator<<(CGAL::Window_stream& os, const X_curve& c){
 leda_segment s(c.xcoord1D(),c.ycoord1D(),c.xcoord2D(),c.ycoord2D()); 
 return os << s; 
 }
 
  #endif
*/
/*
#define UNBOUNDED_CURVE			100
#define TARGET_UNBOUNDED_CURVE	101
#define SOURCE_UNBOUNDED_CURVE	110
#define BOUNDED_CURVE			111
*/
#define UNBOUNDED_X_CURVE "Line"
#define TARGET_UNBOUNDED_X_CURVE "Ray"
#define BOUNDED_X_CURVE "Segment"

#define INSERT_ACTION 0
#define REMOVE_ACTION 1
#define SPLIT_ACTION 2
/*
#define MERGE_ACTION 3 
#define LOCATE_ACTION 4 
#define VERTICAL_RAY_SHOOT_ACTION 5 
*/
#define LOCATE_ACTION 3 
#define VERTICAL_RAY_SHOOT_ACTION 4 

#define REFRESH_BUTTON 2
#define REMOVE_ALL_BUTTON 1
#define EXIT_BUTTON 0

#define LINE_TYPE 0
#define RAY_TYPE 1
#define SEGMENT_TYPE 2
#define POINT_TYPE 3

#define OUTPUT_OPERATOR_MODE 0
#define OUTPUT_OPERATOR_BBOX_MODE 1
#define WRITE_MODE 2
#define WRITE_BBOX_MODE 3

extern int status;
extern bool redraw_status;
extern X_curve_container curve_list;
extern Point_container point_list;

#endif // CGAL_USE_LEDA

#endif //__DRAW_MAP_H

