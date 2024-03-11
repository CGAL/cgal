/**
 * @file   algo/3d/AbstractVertexSplitter.h
 * @author Gernot Walzl
 * @date   2012-10-17
 */

#ifndef ALGO_3D_ABSTRACTVERTEXSPLITTER_H
#define ALGO_3D_ABSTRACTVERTEXSPLITTER_H

#include "algo/3d/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <string>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class AbstractVertexSplitter {
public:
    virtual ~AbstractVertexSplitter();

    virtual PolyhedronSPtr splitVertex(VertexSPtr vertex) = 0;  // abstract

    static const int ANGLE_VERTEX_SPLITTER = -1;  // does not work
    static const int COMBI_VERTEX_SPLITTER = 1;
    static const int CONVEX_VERTEX_SPLITTER = 2;
    static const int VOLUME_VERTEX_SPLITTER = 3;
    static const int WEIGHT_VERTEX_SPLITTER = 4;
    static const int SPHERE_VERTEX_SPLITTER = 5;

    virtual int getType() const;

    static PolyhedronSPtr splitConvexVertex(VertexSPtr vertex);
    static PolyhedronSPtr splitReflexVertex(VertexSPtr vertex);

    static PolyhedronSPtr shiftFacets(PolyhedronSPtr polyhedron, double offset);
    static bool checkSplitted(PolyhedronSPtr polyhedron);

    virtual std::string toString() const;

protected:
    AbstractVertexSplitter();

    int type_;
};

} }

#endif /* ALGO_3D_ABSTRACTVERTEXSPLITTER_H */
