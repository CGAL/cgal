/**
 * @file   algo/3d/ConvexVertexSplitter.h
 * @author Gernot Walzl
 * @date   2013-02-01
 */

#ifndef ALGO_3D_CONVEXVERTEXSPLITTER_H
#define ALGO_3D_CONVEXVERTEXSPLITTER_H

#include "algo/ptrs.h"
#include "algo/3d/ptrs.h"
#include "algo/3d/CombiVertexSplitter.h"
#include "data/3d/ptrs.h"
#include <string>

namespace algo { namespace _3d {

class ConvexVertexSplitter : public CombiVertexSplitter {
public:
    virtual ~ConvexVertexSplitter();

    static ConvexVertexSplitterSPtr create();

    static int countConvexEdges(PolyhedronSPtr polyhedron);

    virtual PolyhedronSPtr splitVertex(VertexSPtr vertex);

    virtual std::string toString() const;

protected:
    ConvexVertexSplitter();
    int optimization_;
};

} }

#endif /* ALGO_3D_CONVEXVERTEXSPLITTER_H */

