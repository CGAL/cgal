/**
 * @file   data/3d/CircularVertexData.h
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#ifndef DATA_3D_CIRCULARVERTEXDATA_H
#define DATA_3D_CIRCULARVERTEXDATA_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

/*!
 * This class is intended to be subclassed.
 */
class CircularVertexData {
public:
    virtual ~CircularVertexData();

    CircularVertexSPtr getVertex() const;
    void setVertex(CircularVertexSPtr vertex);

    bool isHighlight() const;
    void setHighlight(bool highlight);

protected:
    CircularVertexData();
    CircularVertexWPtr vertex_;
    bool highlight_;
};

} }

#endif /* DATA_3D_CIRCULARVERTEXDATA_H */

