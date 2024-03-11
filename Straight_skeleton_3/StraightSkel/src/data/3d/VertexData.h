/**
 * @file   data/3d/VertexData.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_3D_VERTEXDATA_H
#define DATA_3D_VERTEXDATA_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

/*!
 * This class is intended to be subclassed.
 */
class VertexData {
public:
    virtual ~VertexData();

    static VertexDataSPtr create(VertexSPtr vertex);

    VertexSPtr getVertex() const;
    void setVertex(VertexSPtr vertex);

    bool isHighlight() const;
    void setHighlight(bool highlight);

protected:
    VertexData();
    VertexWPtr vertex_;
    bool highlight_;
};

} }

#endif /* DATA_3D_VERTEXDATA_H */
