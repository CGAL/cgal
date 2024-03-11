/**
 * @file   data/2d/skel/ptrs.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_2D_SKEL_PTRS_H
#define DATA_2D_SKEL_PTRS_H

#include "smarter_ptr.h"

/*
 * forward declare classes
 * to break circular includes
 *
 * http://www.boost.org/doc/libs/release/libs/smart_ptr/smart_ptr.htm
 */
namespace data { namespace _2d { namespace skel {

class StraightSkeleton;
class AbstractEvent;
class Node;
class Arc;

class ConstOffsetEvent;
class EdgeEvent;
class SplitEvent;
class TriangleEvent;

class SkelVertexData;
class SkelEdgeData;

typedef SHARED_PTR<StraightSkeleton> StraightSkeletonSPtr;
typedef WEAK_PTR<StraightSkeleton> StraightSkeletonWPtr;
typedef SHARED_PTR<AbstractEvent> AbstractEventSPtr;
typedef WEAK_PTR<AbstractEvent> AbstractEventWPtr;
typedef SHARED_PTR<Node> NodeSPtr;
typedef WEAK_PTR<Node> NodeWPtr;
typedef SHARED_PTR<Arc> ArcSPtr;
typedef WEAK_PTR<Arc> ArcWPtr;

typedef SHARED_PTR<ConstOffsetEvent> ConstOffsetEventSPtr;
typedef WEAK_PTR<ConstOffsetEvent> ConstOffsetEventWPtr;
typedef SHARED_PTR<EdgeEvent> EdgeEventSPtr;
typedef WEAK_PTR<EdgeEvent> EdgeEventWPtr;
typedef SHARED_PTR<SplitEvent> SplitEventSPtr;
typedef WEAK_PTR<SplitEvent> SplitEventWPtr;
typedef SHARED_PTR<TriangleEvent> TriangleEventSPtr;
typedef WEAK_PTR<TriangleEvent> TriangleEventWPtr;

typedef SHARED_PTR<SkelVertexData> SkelVertexDataSPtr;
typedef WEAK_PTR<SkelVertexData> SkelVertexDataWPtr;
typedef SHARED_PTR<SkelEdgeData> SkelEdgeDataSPtr;
typedef WEAK_PTR<SkelEdgeData> SkelEdgeDataWPtr;

} } }

#endif /* DATA_2D_SKEL_PTRS_H */
