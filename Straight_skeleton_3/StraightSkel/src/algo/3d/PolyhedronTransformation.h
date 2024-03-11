/**
 * @file   algo/3d/PolyhedronTransformation.h
 * @author Gernot Walzl
 * @date   2012-09-01
 */

#ifndef ALGO_3D_POLYHEDRONTRANSFORMATION_H
#define ALGO_3D_POLYHEDRONTRANSFORMATION_H

#include "data/3d/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class PolyhedronTransformation {
public:
    virtual ~PolyhedronTransformation();

    static void translate(PolyhedronSPtr polyhedron, Vector3SPtr t);
    static void scale(PolyhedronSPtr polyhedron, Vector3SPtr s);

    /**
     * To check for parallel planes is not enough.
     */
    static bool hasParallelPlanes(PolyhedronSPtr polyhedron);
    static bool doAll3PlanesIntersect(PolyhedronSPtr polyhedron);
    static void randMovePoints(PolyhedronSPtr polyhedron, double range);

    static Point3SPtr boundingBoxMin(PolyhedronSPtr polyhedron);
    static Point3SPtr boundingBoxMax(PolyhedronSPtr polyhedron);

    static void translateNscale(PolyhedronSPtr polyhedron,
            Point3SPtr p_box_min, Point3SPtr p_box_max);

    static bool isInsideBox(PolyhedronSPtr polyhedron,
            Point3SPtr p_box_min, Point3SPtr p_box_max);

protected:
    static Vector3SPtr randVec(double min, double max);
    PolyhedronTransformation();
};

} }

#endif /* ALGO_3D_POLYHEDRONTRANSFORMATION_H */
