#ifndef _MK___TYPEDEFS_H
#define _MK___TYPEDEFS_H

#include <CGAL/basic.h>

//#include <CGAL/Filtered_exact.h>
//typedef CGAL::Filtered_exact<double,CGAL::Gmpq> NT;

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>

#include <CGAL/Segment_Voronoi_diagram_2.h>
#include <CGAL/Segment_Voronoi_diagram_hierarchy_2.h>
#include <CGAL/Segment_Voronoi_diagram_filtered_traits_2.h>

struct Rep : public CGAL::Simple_cartesian<double> {};
//struct Rep : public CGAL::Simple_cartesian<CGAL::Gmpq> {};
//struct Rep : public CGAL::Simple_cartesian<NT> {};

namespace CGAL {
  CGAL::Gmpq sqrt(const CGAL::Gmpq& x)
  {
    return CGAL::Gmpq( sqrt(x.to_double()) );
  }
}

struct Gt
  : public CGAL::Segment_Voronoi_diagram_filtered_traits_2<Rep> {};

//struct Gt
//  : public CGAL::Segment_Voronoi_diagram_traits_2<Rep,CGAL::Ring_tag> {};


typedef Gt::Point_2            Point;
typedef Gt::Segment_2          Segment;
typedef CGAL::Polygon_2<Rep>   Polygon_2;
typedef Gt::Site_2             Site;

//typedef CGAL::Segment_Voronoi_diagram_hierarchy_2<Gt>          SVD_2;
typedef CGAL::Segment_Voronoi_diagram_2<Gt>          SVD_2;

#endif  // _MK___TYPEDEFS_H
