/**
 * @file   algo/3d/ptrs.h
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#ifndef ALGO_3D_PTRS_H
#define ALGO_3D_PTRS_H

#include "smarter_ptr.h"

namespace algo { namespace _3d {

class SimpleStraightSkel;
class AbstractVertexSplitter;
class AngleVertexSplitter;
class CombiVertexSplitter;
class ConvexVertexSplitter;
class VolumeVertexSplitter;
class WeightVertexSplitter;
class SphereVertexSplitter;

class AbstractSimpleSphericalSkel;
class ProjSimpleSphericalSkel;
class RotSimpleSphericalSkel;
class TransSimpleSphericalSkel;
class SpeedSimpleSphericalSkel;

class GraphChecker;

typedef SHARED_PTR<SimpleStraightSkel> SimpleStraightSkelSPtr;
typedef WEAK_PTR<SimpleStraightSkel> SimpleStraightSkelWPtr;
typedef SHARED_PTR<AbstractVertexSplitter> AbstractVertexSplitterSPtr;
typedef WEAK_PTR<AbstractVertexSplitter> AbstractVertexSplitterWPtr;
typedef SHARED_PTR<AngleVertexSplitter> AngleVertexSplitterSPtr;
typedef WEAK_PTR<AngleVertexSplitter> AngleVertexSplitterWPtr;
typedef SHARED_PTR<CombiVertexSplitter> CombiVertexSplitterSPtr;
typedef WEAK_PTR<CombiVertexSplitter> CombiVertexSplitterWPtr;
typedef SHARED_PTR<ConvexVertexSplitter> ConvexVertexSplitterSPtr;
typedef WEAK_PTR<ConvexVertexSplitter> ConvexVertexSplitterWPtr;
typedef SHARED_PTR<VolumeVertexSplitter> VolumeVertexSplitterSPtr;
typedef WEAK_PTR<VolumeVertexSplitter> VolumeVertexSplitterWPtr;
typedef SHARED_PTR<WeightVertexSplitter> WeightVertexSplitterSPtr;
typedef WEAK_PTR<WeightVertexSplitter> WeightVertexSplitterWPtr;
typedef SHARED_PTR<SphereVertexSplitter> SphereVertexSplitterSPtr;
typedef WEAK_PTR<SphereVertexSplitter> SphereVertexSplitterWPtr;

typedef SHARED_PTR<AbstractSimpleSphericalSkel> AbstractSimpleSphericalSkelSPtr;
typedef WEAK_PTR<AbstractSimpleSphericalSkel> AbstractSimpleSphericalSkelWPtr;
typedef SHARED_PTR<ProjSimpleSphericalSkel> ProjSimpleSphericalSkelSPtr;
typedef WEAK_PTR<ProjSimpleSphericalSkel> ProjSimpleSphericalSkelWPtr;
typedef SHARED_PTR<RotSimpleSphericalSkel> RotSimpleSphericalSkelSPtr;
typedef WEAK_PTR<RotSimpleSphericalSkel> RotSimpleSphericalSkelWPtr;
typedef SHARED_PTR<TransSimpleSphericalSkel> TransSimpleSphericalSkelSPtr;
typedef WEAK_PTR<TransSimpleSphericalSkel> TransSimpleSphericalSkelWPtr;
typedef SHARED_PTR<SpeedSimpleSphericalSkel> SpeedSimpleSphericalSkelSPtr;
typedef WEAK_PTR<SpeedSimpleSphericalSkel> SpeedSimpleSphericalSkelWPtr;

typedef SHARED_PTR<GraphChecker> GraphCheckerSPtr;
typedef WEAK_PTR<GraphChecker> GraphCheckerWPtr;

} }

#endif /* ALGO_3D_PTRS_H */
