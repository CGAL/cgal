/**
 * @file   data/3d/skel/ptrs.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_SKEL_PTRS_H
#define DATA_3D_SKEL_PTRS_H

#include "smarter_ptr.h"

/*
 * forward declare classes
 * to break circular includes
 *
 * http://www.boost.org/doc/libs/release/libs/smart_ptr/smart_ptr.htm
 */
namespace data { namespace _3d { namespace skel {

class StraightSkeleton;
class Node;
class Arc;
class Sheet;

class AbstractEvent;
class ConstOffsetEvent;
class SaveOffsetEvent;
class EdgeEvent;
class EdgeMergeEvent;
class TriangleEvent;
class DblEdgeMergeEvent;
class DblTriangleEvent;
class TetrahedronEvent;
class VertexEvent;
class FlipVertexEvent;
class SurfaceEvent;
class PolyhedronSplitEvent;
class SplitMergeEvent;
class EdgeSplitEvent;
class PierceEvent;

class SkelEdgeData;
class SkelVertexData;
class SkelFacetData;

class SphericalSkeleton;
class CircularNode;
class CircularArc;

class SphericalAbstractEvent;
class SphericalConstOffsetEvent;
class SphericalEdgeEvent;
class SphericalSplitEvent;
class SphericalTriangleEvent;
class SphericalDblEdgeEvent;
class SphericalLeaveEvent;
class SphericalReturnEvent;
class SphericalDblLeaveEvent;
class SphericalDblReturnEvent;
class SphericalVertexEvent;
class SphericalEdgeMergeEvent;
class SphericalInversionEvent;

class SphericalSkelEdgeData;
class SphericalSkelVertexData;

class SphericalOffset;

typedef SHARED_PTR<StraightSkeleton> StraightSkeletonSPtr;
typedef WEAK_PTR<StraightSkeleton> StraightSkeletonWPtr;
typedef SHARED_PTR<Node> NodeSPtr;
typedef WEAK_PTR<Node> NodeWPtr;
typedef SHARED_PTR<Arc> ArcSPtr;
typedef WEAK_PTR<Arc> ArcWPtr;
typedef SHARED_PTR<Sheet> SheetSPtr;
typedef WEAK_PTR<Sheet> SheetWPtr;

typedef SHARED_PTR<AbstractEvent> AbstractEventSPtr;
typedef WEAK_PTR<AbstractEvent> AbstractEventWPtr;
typedef SHARED_PTR<ConstOffsetEvent> ConstOffsetEventSPtr;
typedef WEAK_PTR<ConstOffsetEvent> ConstOffsetEventWPtr;
typedef SHARED_PTR<SaveOffsetEvent> SaveOffsetEventSPtr;
typedef WEAK_PTR<SaveOffsetEvent> SaveOffsetEventWPtr;
typedef SHARED_PTR<EdgeEvent> EdgeEventSPtr;
typedef WEAK_PTR<EdgeEvent> EdgeEventWPtr;
typedef SHARED_PTR<EdgeMergeEvent> EdgeMergeEventSPtr;
typedef WEAK_PTR<EdgeMergeEvent> EdgeMergeEventWPtr;
typedef SHARED_PTR<TriangleEvent> TriangleEventSPtr;
typedef WEAK_PTR<TriangleEvent> TriangleEventWPtr;
typedef SHARED_PTR<DblEdgeMergeEvent> DblEdgeMergeEventSPtr;
typedef WEAK_PTR<DblEdgeMergeEvent> DblEdgeMergeEventWPtr;
typedef SHARED_PTR<DblTriangleEvent> DblTriangleEventSPtr;
typedef WEAK_PTR<DblTriangleEvent> DblTriangleEventWPtr;
typedef SHARED_PTR<TetrahedronEvent> TetrahedronEventSPtr;
typedef WEAK_PTR<TetrahedronEvent> TetrahedronEventWPtr;
typedef SHARED_PTR<VertexEvent> VertexEventSPtr;
typedef WEAK_PTR<VertexEvent> VertexEventWPtr;
typedef SHARED_PTR<FlipVertexEvent> FlipVertexEventSPtr;
typedef WEAK_PTR<FlipVertexEvent> FlipVertexEventWPtr;
typedef SHARED_PTR<SurfaceEvent> SurfaceEventSPtr;
typedef WEAK_PTR<SurfaceEvent> SurfaceEventWPtr;
typedef SHARED_PTR<PolyhedronSplitEvent> PolyhedronSplitEventSPtr;
typedef WEAK_PTR<PolyhedronSplitEvent> PolyhedronSplitEventWPtr;
typedef SHARED_PTR<SplitMergeEvent> SplitMergeEventSPtr;
typedef WEAK_PTR<SplitMergeEvent> SplitMergeEventWPtr;
typedef SHARED_PTR<EdgeSplitEvent> EdgeSplitEventSPtr;
typedef WEAK_PTR<EdgeSplitEvent> EdgeSplitEventWPtr;
typedef SHARED_PTR<PierceEvent> PierceEventSPtr;
typedef WEAK_PTR<PierceEvent> PierceEventWPtr;

typedef SHARED_PTR<SkelEdgeData> SkelEdgeDataSPtr;
typedef WEAK_PTR<SkelEdgeData> SkelEdgeDataWPtr;
typedef SHARED_PTR<SkelVertexData> SkelVertexDataSPtr;
typedef WEAK_PTR<SkelVertexData> SkelVertexDataWPtr;
typedef SHARED_PTR<SkelFacetData> SkelFacetDataSPtr;
typedef WEAK_PTR<SkelFacetData> SkelFacetDataWPtr;

typedef SHARED_PTR<SphericalSkeleton> SphericalSkeletonSPtr;
typedef WEAK_PTR<SphericalSkeleton> SphericalSkeletonWPtr;
typedef SHARED_PTR<CircularNode> CircularNodeSPtr;
typedef WEAK_PTR<CircularNode> CircularNodeWPtr;
typedef SHARED_PTR<CircularArc> CircularArcSPtr;
typedef WEAK_PTR<CircularArc> CircularArcWPtr;

typedef SHARED_PTR<SphericalAbstractEvent> SphericalAbstractEventSPtr;
typedef WEAK_PTR<SphericalAbstractEvent> SphericalAbstractEventWPtr;
typedef SHARED_PTR<SphericalConstOffsetEvent> SphericalConstOffsetEventSPtr;
typedef WEAK_PTR<SphericalConstOffsetEvent> SphericalConstOffsetEventWPtr;
typedef SHARED_PTR<SphericalEdgeEvent> SphericalEdgeEventSPtr;
typedef WEAK_PTR<SphericalEdgeEvent> SphericalEdgeEventWPtr;
typedef SHARED_PTR<SphericalSplitEvent> SphericalSplitEventSPtr;
typedef WEAK_PTR<SphericalSplitEvent> SphericalSplitEventWPtr;
typedef SHARED_PTR<SphericalTriangleEvent> SphericalTriangleEventSPtr;
typedef WEAK_PTR<SphericalTriangleEvent> SphericalTriangleEventWPtr;
typedef SHARED_PTR<SphericalDblEdgeEvent> SphericalDblEdgeEventSPtr;
typedef WEAK_PTR<SphericalDblEdgeEvent> SphericalDblEdgeEventWPtr;
typedef SHARED_PTR<SphericalLeaveEvent> SphericalLeaveEventSPtr;
typedef WEAK_PTR<SphericalLeaveEvent> SphericalLeaveEventWPtr;
typedef SHARED_PTR<SphericalReturnEvent> SphericalReturnEventSPtr;
typedef WEAK_PTR<SphericalReturnEvent> SphericalReturnEventWPtr;
typedef SHARED_PTR<SphericalDblLeaveEvent> SphericalDblLeaveEventSPtr;
typedef WEAK_PTR<SphericalDblLeaveEvent> SphericalDblLeaveEventWPtr;
typedef SHARED_PTR<SphericalDblReturnEvent> SphericalDblReturnEventSPtr;
typedef WEAK_PTR<SphericalDblReturnEvent> SphericalDblReturnEventWPtr;
typedef SHARED_PTR<SphericalVertexEvent> SphericalVertexEventSPtr;
typedef WEAK_PTR<SphericalVertexEvent> SphericalVertexEventWPtr;
typedef SHARED_PTR<SphericalEdgeMergeEvent> SphericalEdgeMergeEventSPtr;
typedef WEAK_PTR<SphericalEdgeMergeEvent> SphericalEdgeMergeEventWPtr;
typedef SHARED_PTR<SphericalInversionEvent> SphericalInversionEventSPtr;
typedef WEAK_PTR<SphericalInversionEvent> SphericalInversionEventWPtr;

typedef SHARED_PTR<SphericalSkelEdgeData> SphericalSkelEdgeDataSPtr;
typedef WEAK_PTR<SphericalSkelEdgeData> SphericalSkelEdgeDataWPtr;
typedef SHARED_PTR<SphericalSkelVertexData> SphericalSkelVertexDataSPtr;
typedef WEAK_PTR<SphericalSkelVertexData> SphericalSkelVertexDataWPtr;

typedef SHARED_PTR<SphericalOffset> SphericalOffsetSPtr;
typedef WEAK_PTR<SphericalOffset> SphericalOffsetWPtr;

} } }

#endif /* DATA_3D_SKEL_PTRS_H */
