#ifndef __DEMO_H
#define __DEMO_H

// if LEDA is not installed, a message will be issued in runtime by demo.C.
#ifdef CGAL_USE_LEDA

#include "configuration"

/**************************************
Configuring bounding box 
***************************************/

#if defined(USE_DYNAMIC_CLOSED_BOUNDING_BOX)

#define BOUNDING_BOX CGAL::Pm_dynamic_closed_bounding_box<Planar_map>
#include <CGAL/Pm_dynamic_closed_bounding_box.h>

#else 
#if defined(USE_DYNAMIC_OPEN_BOUNDING_BOX)

#define BOUNDING_BOX CGAL::Pm_dynamic_open_bounding_box<Planar_map>
#include <CGAL/Pm_dynamic_open_bounding_box.h>

#else
#if defined(USE_STATIC_CLOSED_BOUNDING_BOX)

#define BOUNDING_BOX CGAL::Pm_static_closed_bounding_box<Planar_map>
#include <CGAL/Pm_static_closed_bounding_box.h>

#else
#if defined(USE_STATIC_OPEN_BOUNDING_BOX)

#define BOUNDING_BOX CGAL::Pm_static_open_bounding_box<Planar_map>
#include <CGAL/Pm_static_open_bounding_box.h>

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

typedef CGAL::Window_stream  Window_stream;

#endif // CGAL_USE_LEDA

#endif // __DEMO_H

