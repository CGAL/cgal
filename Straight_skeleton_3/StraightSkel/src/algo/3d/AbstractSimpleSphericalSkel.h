/**
 * @file   algo/3d/AbstractSimpleSphericalSkel.h
 * @author Gernot Walzl
 * @date   2013-01-11
 */

#ifndef ALGO_3D_ABSTRACTSIMPLESPHERICALSKEL_H
#define ALGO_3D_ABSTRACTSIMPLESPHERICALSKEL_H

#include "typedefs_thread.h"
#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class AbstractSimpleSphericalSkel {
public:
    virtual ~AbstractSimpleSphericalSkel();

    static const int PROJ_SIMPLE_SPHERICAL_SKEL = 1;
    static const int ROT_SIMPLE_SPHERICAL_SKEL = 2;
    static const int TRANS_SIMPLE_SPHERICAL_SKEL = 3;
    static const int SPEED_SIMPLE_SPHERICAL_SKEL = 4;

    virtual int getType() const;

    virtual void run() = 0;  // abstract
    virtual ThreadSPtr startThread();

    virtual bool isReflex(CircularVertexSPtr vertex) = 0;  // abstract

    virtual bool init(SphericalPolygonSPtr polygon) = 0;  // abstract

    static CircularEdgeSPtr getEdgeOrigin(CircularEdgeSPtr edge);
    virtual CircularArcSPtr createArc(CircularVertexSPtr vertex);

    static Point3SPtr vanishesAt(CircularEdgeSPtr edge);
    static Point3SPtr crashAt(CircularVertexSPtr vertex, CircularEdgeSPtr edge);

    static bool isTriangle(CircularEdgeSPtr edge_begin);

    virtual void appendEventNode(CircularNodeSPtr node);

    virtual SphericalSkeletonSPtr getResult() const;

protected:
    AbstractSimpleSphericalSkel();

    int type_;
    SphericalPolygonSPtr polygon_;
    ControllerSPtr controller_;
    SphericalSkeletonSPtr skel_result_;
};

} }

#endif /* ALGO_3D_ABSTRACTSIMPLESPHERICALSKEL_H */
