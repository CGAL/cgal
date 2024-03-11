/**
 * @file   algo/3d/SphereVertexSplitter.h
 * @author Gernot Walzl
 * @date   2012-11-20
 */

#ifndef ALGO_3D_SPHEREVERTEXSPLITTER_H
#define ALGO_3D_SPHEREVERTEXSPLITTER_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractVertexSplitter.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalSkeleton.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using data::_3d::skel::SphericalSkeleton;
using data::_3d::skel::SphericalSkeletonSPtr;

class SphereVertexSplitter : public AbstractVertexSplitter {
public:
    virtual ~SphereVertexSplitter();

    static SphereVertexSplitterSPtr create();
    static SphereVertexSplitterSPtr create(ControllerSPtr controller);

    static FacetSPtr findFacet(EdgeSPtr edge_1, EdgeSPtr edge_2);

    static SphericalPolygonSPtr intersectPolyhedron(VertexSPtr vertex);

    virtual PolyhedronSPtr splitVertex(VertexSPtr vertex);

protected:
    SphereVertexSplitter();
    SphereVertexSplitter(ControllerSPtr controller);
    ControllerSPtr controller_;
};

} }

#endif /* ALGO_3D_SPHEREVERTEXSPLITTER_H */
