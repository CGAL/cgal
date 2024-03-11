/**
 * @file   data/2d/ptrs.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_2D_PTRS_H
#define DATA_2D_PTRS_H

#include "smarter_ptr.h"

#include "config.h"
#ifdef USE_CGAL
    #include "cgal_kernel.h"
#else
    #include "kernel/Point2.h"
    #include "kernel/Line2.h"
    #include "kernel/Segment2.h"
    #include "kernel/Vector2.h"
#endif

/*
 * forward declare classes
 * to break circular includes
 *
 * http://www.boost.org/doc/libs/release/libs/smart_ptr/smart_ptr.htm
 */
namespace data { namespace _2d {

#ifdef USE_CGAL
using CGAL::Point2;
using CGAL::Line2;
using CGAL::Segment2;
using CGAL::Vector2;
#else
using kernel::Point2;
using kernel::Line2;
using kernel::Segment2;
using kernel::Vector2;
#endif

class Polygon;
class Vertex;
class Edge;
class VertexData;
class EdgeData;

typedef SHARED_PTR<Point2> Point2SPtr;
typedef WEAK_PTR<Point2> Point2WPtr;
typedef SHARED_PTR<Line2> Line2SPtr;
typedef WEAK_PTR<Line2> Line2WPtr;
typedef SHARED_PTR<Segment2> Segment2SPtr;
typedef WEAK_PTR<Segment2> Segment2WPtr;
typedef SHARED_PTR<Vector2> Vector2SPtr;
typedef WEAK_PTR<Vector2> Vector2WPtr;

typedef SHARED_PTR<Polygon> PolygonSPtr;
typedef WEAK_PTR<Polygon> PolygonWPtr;
typedef SHARED_PTR<Vertex> VertexSPtr;
typedef WEAK_PTR<Vertex> VertexWPtr;
typedef SHARED_PTR<Edge> EdgeSPtr;
typedef WEAK_PTR<Edge> EdgeWPtr;
typedef SHARED_PTR<VertexData> VertexDataSPtr;
typedef WEAK_PTR<VertexData> VertexDataWPtr;
typedef SHARED_PTR<EdgeData> EdgeDataSPtr;
typedef WEAK_PTR<EdgeData> EdgeDataWPtr;

} }

#endif /* DATA_2D_PTRS_H */
