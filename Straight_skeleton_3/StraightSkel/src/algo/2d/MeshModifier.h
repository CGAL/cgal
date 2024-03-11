/**
 * @file   algo/2d/MeshModifier.h
 * @author Gernot Walzl
 * @date   2014-04-27
 */

#ifndef ALGO_2D_MESHMODIFIER_H
#define ALGO_2D_MESHMODIFIER_H

#include "data/2d/ptrs.h"
#include "data/2d/mesh/ptrs.h"

namespace algo { namespace _2d {

using namespace data::_2d;
using namespace data::_2d::mesh;

class MeshModifier {
public:
    virtual ~MeshModifier();

    static void splitEdge(MeshVertexSPtr v_src, MeshVertexSPtr v_dst,
            MeshVertexSPtr v_insert);
    static MeshCellSPtr splitCell(MeshCellSPtr cell,
            MeshVertexSPtr v_src, MeshVertexSPtr v_dst);
    static void mergeCells(MeshCellSPtr cell_1, MeshCellSPtr cell_2);
    static void mergeVertices(MeshVertexSPtr vertex_1, MeshVertexSPtr vertex_2);

protected:
    MeshModifier();
};

} }

#endif /* ALGO_2D_MESHMODIFIER_H */
