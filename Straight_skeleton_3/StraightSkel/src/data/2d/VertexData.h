/**
 * @file   data/2d/VertexData.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_2D_VERTEXDATA_H
#define DATA_2D_VERTEXDATA_H

#include "data/2d/ptrs.h"

namespace data { namespace _2d {

/*!
 * This class is intended to be subclassed.
 */
class VertexData {
public:
    virtual ~VertexData();

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

#endif /* DATA_2D_VERTEXDATA_H */

