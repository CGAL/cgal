#ifndef _MK___TYPEDEFS_H
#define _MK___TYPEDEFS_H

#include <CGAL/basic.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_traits_2.h>

struct Rep : public CGAL::Simple_cartesian<double> {};

struct Gt
  : public CGAL::Segment_Voronoi_diagram_filtered_traits_2<Rep> {};


typedef Gt::Point_2            Point;
typedef Gt::Segment_2          Segment;
typedef CGAL::Polygon_2<Rep>   Polygon;
typedef Gt::Site_2             Site;

typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>          SVD_2;

#endif  // _MK___TYPEDEFS_H
