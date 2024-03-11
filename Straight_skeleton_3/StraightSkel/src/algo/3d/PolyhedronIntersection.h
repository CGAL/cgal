/**
 * @file   algo/3d/PolyhedronIntersection.h
 * @author Gernot Walzl
 * @date   2012-09-26
 */

#ifndef ALGO_3D_POLYHEDRONINTERSECTION_H
#define ALGO_3D_POLYHEDRONINTERSECTION_H

#include "data/2d/ptrs.h"
#include "data/3d/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class PolyhedronIntersection {
public:
    virtual ~PolyhedronIntersection();

    static Point3SPtr intersect(EdgeSPtr edge, Plane3SPtr plane);
    static EdgeSPtr findDst(FacetSPtr facet, EdgeSPtr edge_src,
            Plane3SPtr plane);
    static FacetSPtr intersect(PolyhedronSPtr polyhedron, Plane3SPtr plane);

    static data::_2d::PolygonSPtr toWeighted2d(FacetSPtr facet);

protected:
    PolyhedronIntersection();
};

} }

#endif /* ALGO_3D_POLYHEDRONINTERSECTION_H */
