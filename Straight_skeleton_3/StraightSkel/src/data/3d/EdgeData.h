/**
 * @file   data/3d/EdgeData.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_3D_EDGEDATA_H
#define DATA_3D_EDGEDATA_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

/*!
 * This class is intended to be subclassed.
 */
class EdgeData {
public:
    virtual ~EdgeData();

    static EdgeDataSPtr create(EdgeSPtr edge);

    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);

    bool isHighlight() const;
    void setHighlight(bool highlight);

protected:
    EdgeData();
    EdgeWPtr edge_;
    bool highlight_;
};

} }

#endif /* DATA_3D_EDGEDATA_H */

