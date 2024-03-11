/**
 * @file   data/3d/ptrs.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_PTRS_H
#define DATA_3D_PTRS_H

#include "smarter_ptr.h"

#include "config.h"
#ifdef USE_CGAL
    #include "cgal_kernel.h"
#else
    #include "kernel/Point3.h"
    #include "kernel/Vector3.h"
    #include "kernel/Line3.h"
    #include "kernel/Segment3.h"
    #include "kernel/Plane3.h"
    #include "kernel/Sphere3.h"
#endif

/*
 * forward declare classes
 * to break circular includes
 *
 * http://www.boost.org/doc/libs/release/libs/smart_ptr/smart_ptr.htm
 */
namespace data { namespace _3d {

#ifdef USE_CGAL
using CGAL::Point3;
using CGAL::Vector3;
using CGAL::Line3;
using CGAL::Segment3;
using CGAL::Plane3;
using CGAL::Sphere3;
#else
using kernel::Point3;
using kernel::Vector3;
using kernel::Line3;
using kernel::Segment3;
using kernel::Plane3;
using kernel::Sphere3;
#endif

class Polyhedron;
class Facet;
class FacetData;
class Vertex;
class VertexData;
class Edge;
class EdgeData;
class Triangle;

class SphericalPolygon;
class CircularVertex;
class CircularVertexData;
class CircularEdge;
class CircularEdgeData;

typedef SHARED_PTR<Point3> Point3SPtr;
typedef WEAK_PTR<Point3> Point3WPtr;
typedef SHARED_PTR<Vector3> Vector3SPtr;
typedef WEAK_PTR<Vector3> Vector3WPtr;
typedef SHARED_PTR<Line3> Line3SPtr;
typedef WEAK_PTR<Line3> Line3WPtr;
typedef SHARED_PTR<Segment3> Segment3SPtr;
typedef WEAK_PTR<Segment3> Segment3WPtr;
typedef SHARED_PTR<Plane3> Plane3SPtr;
typedef WEAK_PTR<Plane3> Plane3WPtr;
typedef SHARED_PTR<Sphere3> Sphere3SPtr;
typedef WEAK_PTR<Sphere3> Sphere3WPtr;

typedef SHARED_PTR<Polyhedron> PolyhedronSPtr;
typedef WEAK_PTR<Polyhedron> PolyhedronWPtr;
typedef SHARED_PTR<Facet> FacetSPtr;
typedef WEAK_PTR<Facet> FacetWPtr;
typedef SHARED_PTR<FacetData> FacetDataSPtr;
typedef WEAK_PTR<FacetData> FacetDataWPtr;
typedef SHARED_PTR<Vertex> VertexSPtr;
typedef WEAK_PTR<Vertex> VertexWPtr;
typedef SHARED_PTR<VertexData> VertexDataSPtr;
typedef WEAK_PTR<VertexData> VertexDataWPtr;
typedef SHARED_PTR<Edge> EdgeSPtr;
typedef WEAK_PTR<Edge> EdgeWPtr;
typedef SHARED_PTR<EdgeData> EdgeDataSPtr;
typedef WEAK_PTR<EdgeData> EdgeDataWPtr;
typedef SHARED_PTR<Triangle> TriangleSPtr;
typedef WEAK_PTR<Triangle> TriangleWPtr;

typedef SHARED_PTR<SphericalPolygon> SphericalPolygonSPtr;
typedef WEAK_PTR<SphericalPolygon> SphericalPolygonWPtr;
typedef SHARED_PTR<CircularVertex> CircularVertexSPtr;
typedef WEAK_PTR<CircularVertex> CircularVertexWPtr;
typedef SHARED_PTR<CircularVertexData> CircularVertexDataSPtr;
typedef WEAK_PTR<CircularVertexData> CircularVertexDataWPtr;
typedef SHARED_PTR<CircularEdge> CircularEdgeSPtr;
typedef WEAK_PTR<CircularEdge> CircularEdgeWPtr;
typedef SHARED_PTR<CircularEdgeData> CircularEdgeDataSPtr;
typedef WEAK_PTR<CircularEdgeData> CircularEdgeDataWPtr;

} }

#endif /* DATA_3D_PTRS_H */
