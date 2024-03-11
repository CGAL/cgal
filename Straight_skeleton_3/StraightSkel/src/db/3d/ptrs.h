#ifndef DB_3D_PTRS_H
#define DB_3D_PTRS_H

#include "smarter_ptr.h"

namespace db { namespace _3d {

class PointDAO;
class PlaneDAO;

class VertexDAO;
class EdgeDAO;
class TriangleDAO;
class FacetDAO;
class PolyhedronDAO;

class NodeDAO;
class ArcDAO;
class SheetDAO;
class EventDAO;
class StraightSkeletonDAO;

typedef SHARED_PTR<PointDAO> PointDAOSPtr;
typedef WEAK_PTR<PointDAO> PointDAOWPtr;
typedef SHARED_PTR<PlaneDAO> PlaneDAOSPtr;
typedef WEAK_PTR<PlaneDAO> PlaneDAOWPtr;

typedef SHARED_PTR<VertexDAO> VertexDAOSPtr;
typedef WEAK_PTR<VertexDAO> VertexDAOWPtr;
typedef SHARED_PTR<EdgeDAO> EdgeDAOSPtr;
typedef WEAK_PTR<EdgeDAO> EdgeDAOWPtr;
typedef SHARED_PTR<TriangleDAO> TriangleDAOSPtr;
typedef WEAK_PTR<TriangleDAO> TriangleDAOWPtr;
typedef SHARED_PTR<FacetDAO> FacetDAOSPtr;
typedef WEAK_PTR<FacetDAO> FacetDAOWPtr;
typedef SHARED_PTR<PolyhedronDAO> PolyhedronDAOSPtr;
typedef WEAK_PTR<PolyhedronDAO> PolyhedronDAOWPtr;

typedef SHARED_PTR<NodeDAO> NodeDAOSPtr;
typedef WEAK_PTR<NodeDAO> NodeDAOWPtr;
typedef SHARED_PTR<ArcDAO> ArcDAOSPtr;
typedef WEAK_PTR<ArcDAO> ArcDAOWPtr;
typedef SHARED_PTR<SheetDAO> SheetDAOSPtr;
typedef WEAK_PTR<SheetDAO> SheetDAOWPtr;
typedef SHARED_PTR<EventDAO> EventDAOSPtr;
typedef WEAK_PTR<EventDAO> EventDAOWPtr;
typedef SHARED_PTR<StraightSkeletonDAO> StraightSkeletonDAOSPtr;
typedef WEAK_PTR<StraightSkeletonDAO> StraightSkeletonDAOWPtr;

} }

#endif /* DB_3D_PTRS_H */
