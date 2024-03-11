/**
 * @file   db/2d/ptrs.h
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#ifndef DB_2D_PTRS_H
#define DB_2D_PTRS_H

#include "smarter_ptr.h"

namespace db { namespace _2d {

class PointDAO;

class VertexDAO;
class EdgeDAO;
class PolygonDAO;

class NodeDAO;
class ArcDAO;
class EventDAO;
class StraightSkeletonDAO;

typedef SHARED_PTR<PointDAO> PointDAOSPtr;
typedef WEAK_PTR<PointDAO> PointDAOWPtr;

typedef SHARED_PTR<VertexDAO> VertexDAOSPtr;
typedef WEAK_PTR<VertexDAO> VertexDAOWPtr;
typedef SHARED_PTR<EdgeDAO> EdgeDAOSPtr;
typedef WEAK_PTR<EdgeDAO> EdgeDAOWPtr;
typedef SHARED_PTR<PolygonDAO> PolygonDAOSPtr;
typedef WEAK_PTR<PolygonDAO> PolygonDAOWPtr;

typedef SHARED_PTR<NodeDAO> NodeDAOSPtr;
typedef WEAK_PTR<NodeDAO> NodeDAOWPtr;
typedef SHARED_PTR<ArcDAO> ArcDAOSPtr;
typedef WEAK_PTR<ArcDAO> ArcDAOWPtr;
typedef SHARED_PTR<EventDAO> EventDAOSPtr;
typedef WEAK_PTR<EventDAO> EventDAOWPtr;
typedef SHARED_PTR<StraightSkeletonDAO> StraightSkeletonDAOSPtr;
typedef WEAK_PTR<StraightSkeletonDAO> StraightSkeletonDAOWPtr;

} }

#endif /* DB_2D_PTRS_H */
