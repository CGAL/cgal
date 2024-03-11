/**
 * @file   data/2d/EdgeData.h
 * @author Gernot Walzl
 * @date   2011-11-23
 */

#ifndef DATA_2D_EDGEDATA_H
#define DATA_2D_EDGEDATA_H

#include "data/2d/ptrs.h"

namespace data { namespace _2d {

/*!
 * This class is intended to be subclassed.
 */
class EdgeData {
public:
    virtual ~EdgeData();

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

#endif /* DATA_2D_EDGEDATA_H */
