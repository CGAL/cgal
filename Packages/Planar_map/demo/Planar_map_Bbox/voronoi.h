// if LEDA is not installed, a message will be issued in runtime by demo.C.
#ifdef CGAL_USE_LEDA

#include "configuration"

/**************************************
Configuring bounding box 
***************************************/

#ifdef USE_DYNAMIC_CLOSED_BOUNDING_BOX

#include <CGAL/Pm_dynamic_closed_bounding_box.h>
#define BOUNDING_BOX CGAL::Pm_dynamic_closed_bounding_box<Planar_map>

#else
#ifdef USE_DYNAMIC_OPEN_BOUNDING_BOX

#include <CGAL/Pm_dynamic_open_bounding_box.h>
#define BOUNDING_BOX CGAL::Pm_dynamic_open_bounding_box<Planar_map>

#else 
#ifdef USE_STATIC_CLOSED_BOUNDING_BOX

#include <CGAL/Pm_static_closed_bounding_box.h>
#define BOUNDING_BOX CGAL::Pm_static_closed_bounding_box<Planar_map>

#else
#ifdef USE_STATIC_OPEN_BOUNDING_BOX

#include <CGAL/Pm_static_open_bounding_box.h>
#define BOUNDING_BOX CGAL::Pm_static_open_bounding_box<Planar_map>

#else // USE_UNBOUNDING_BOX

#define BOUNDING_BOX CGAL::Pm_unbounding_box<Planar_map>

#endif
#endif
#endif
#endif

/**************************************
Configuring point location
***************************************/

#ifdef USE_NAIVE_POINT_LOCATION

#include <CGAL/Pm_naive_point_location.h>
#define POINT_LOCATION CGAL::Pm_naive_point_location<Planar_map>
#define POINT_LOCATION_ARGS  
#define POINT_LOCATION_NAME "naive point location"

#else
#if defined(USE_WALK_POINT_LOCATION)

#include <CGAL/Pm_walk_along_line_point_location.h>
#define POINT_LOCATION CGAL::Pm_walk_along_line_point_location<Planar_map>
#define POINT_LOCATION_ARGS  
#define POINT_LOCATION_NAME "walk point location"

#else
#if defined(USE_DEFAULT_WITHOUT_REBUILD)

#define POINT_LOCATION CGAL::Pm_default_point_location<Planar_map>
#define POINT_LOCATION_ARGS (false) 
#define POINT_LOCATION_NAME "default point location (without rebuild)"

#else // USE_DEFAULT_POINT_LOCATION

#define POINT_LOCATION CGAL::Pm_default_point_location<Planar_map>
#define POINT_LOCATION_ARGS (true) 
#define POINT_LOCATION_NAME "default point location"

#endif
#endif
#endif

#include <CGAL/basic.h>
#include <fstream>
#include "draw_map.h"  
#include <CGAL/IO/Planar_map_Window_stream.h>

// Define shorter names to please linker (g++/egcs)
//#define Cartesian Cart

#include <CGAL/Cartesian.h>

#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef CGAL::Window_stream  Window_stream;

#ifdef USE_RATIONAL
typedef leda_rational                          number_type;
#else
typedef double                                 number_type; 
#endif

typedef CGAL::Cartesian<number_type>  Rpst;

typedef CGAL::Point_2<Rpst>  Point;
typedef CGAL::Segment_2<Rpst>  Segment;
typedef CGAL::Ray_2<Rpst>  Ray;
typedef CGAL::Line_2<Rpst>  Line;
typedef CGAL::Triangle_2<Rpst>  Triangle;

typedef CGAL::Triangulation_euclidean_traits_2<Rpst> Gt;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Triangulation;

typedef Triangulation::Face  TFace;
typedef Triangulation::Vertex TVertex;
typedef Triangulation::Face_handle  TFace_handle;
typedef Triangulation::Vertex_handle TVertex_handle;

typedef Triangulation::Face_circulator  TFace_circulator;
typedef Triangulation::Vertex_circulator  TVertex_circulator;

typedef Triangulation::Locate_type TLocate_type;

typedef Triangulation::Face_iterator  TFace_iterator;
typedef Triangulation::Vertex_iterator  TVertex_iterator;
typedef Triangulation::Edge_iterator  TEdge_iterator;
typedef Triangulation::Line_face_circulator  TLine_face_circulator;
/*
typedef std::list<X_curve> X_curve_container;
typedef std::list<Point> Point_container;
*/

#endif // CGAL_USE_LEDA
