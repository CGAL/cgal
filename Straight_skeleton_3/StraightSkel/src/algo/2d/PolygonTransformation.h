/**
 * @file   algo/2d/PolygonTransformation.h
 * @author Gernot Walzl
 * @date   2012-09-18
 */

#ifndef ALGO_2D_POLYGONTRANSFORMATION_H
#define ALGO_2D_POLYGONTRANSFORMATION_H

#include "data/2d/ptrs.h"

namespace algo { namespace _2d {

using namespace data::_2d;

class PolygonTransformation {
public:
    virtual ~PolygonTransformation();

    static bool hasParallelLines(PolygonSPtr polygon);
    static void randMovePoints(PolygonSPtr polygon, double range);

    static Point2SPtr boundingBoxMin(PolygonSPtr polygon);
    static Point2SPtr boundingBoxMax(PolygonSPtr polygon);

protected:
    static Vector2SPtr randVec(double min, double max);
    PolygonTransformation();
};

} }

#endif /* ALGO_2D_POLYGONTRANSFORMATION_H */
