#ifndef BOPS_TYPES_H
#define BOPS_TYPES_H

#include <CGAL/basic.h>
//using std::sqrt;

#include "numrep1.h"
#include <vector>
#include <list>


#include <CGAL/Point_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_2.h>
#include <CGAL/Object.h>

#include "numrep2.h"

using CGAL::Point_2;
using CGAL::Segment_2;
using CGAL::Polygon_2;
using CGAL::Polygon_traits_2;

typedef Point_2<TestRep> TestPoint;
typedef Segment_2<TestRep> TestSegment;

typedef
    Polygon_2<Polygon_traits_2<TestRep>, CGAL_STD::list<TestPoint> >
    TestPolygon;
typedef
    Polygon_2<Polygon_traits_2<TestRep>, CGAL_STD::list<TestPoint> >
    ResultPolygon;

//typedef
//Polygon_2<Polygon_traits_2<TestRep>, CGAL_STD::vector<TestPoint> >
//    TestPolygon;
//typedef
//Polygon_2<Polygon_traits_2<TestRep>, CGAL_STD::vector<TestPoint> >
//    ResultPolygon;

#endif // BOPS_TYPES_H
