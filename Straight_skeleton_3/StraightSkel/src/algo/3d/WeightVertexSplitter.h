/**
 * @file   algo/3d/WeightVertexSplitter.h
 * @author Gernot Walzl
 * @date   2012-10-03
 */

#ifndef ALGO_3D_WEIGHTVERTEXSPLITTER_H
#define ALGO_3D_WEIGHTVERTEXSPLITTER_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/AbstractVertexSplitter.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class WeightVertexSplitter : public AbstractVertexSplitter {
public:
    virtual ~WeightVertexSplitter();

    static WeightVertexSplitterSPtr create();
    static WeightVertexSplitterSPtr create(ControllerSPtr controller);

    static Plane3SPtr createIntersectionPlane(VertexSPtr vertex);
    static FacetSPtr intersectPolyhedron(VertexSPtr vertex, Plane3SPtr plane);
    static data::_2d::PolygonSPtr toWeighted2d(FacetSPtr polygon);

    virtual PolyhedronSPtr splitVertex(VertexSPtr vertex);

protected:
    WeightVertexSplitter();
    WeightVertexSplitter(ControllerSPtr controller);
    ControllerSPtr controller_;
};

} }

#endif /* ALGO_3D_WEIGHTVERTEXSPLITTER_H */
