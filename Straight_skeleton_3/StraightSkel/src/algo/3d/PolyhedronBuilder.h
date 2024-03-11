/**
 * @file   algo/3d/PolyhedronBuilder.h
 * @author Gernot Walzl
 * @date   2012-04-03
 */

#ifndef ALGO_3D_POLYHEDRONBUILDER_H
#define ALGO_3D_POLYHEDRONBUILDER_H

#include "data/3d/ptrs.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class PolyhedronBuilder {
public:
    virtual ~PolyhedronBuilder();
    static PolyhedronSPtr makeTetrahedron(
            Point3SPtr p1, Point3SPtr p2, Point3SPtr p3, Point3SPtr p4);
protected:
    PolyhedronBuilder();
};

} }

#endif /* ALGO_3D_POLYHEDRONBUILDER_H */

