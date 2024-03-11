/**
 * @file   data/3d/CircularEdge.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_CIRCULAREDGE_H
#define DATA_3D_CIRCULAREDGE_H

#include "data/3d/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class CircularEdge {
public:
    virtual ~CircularEdge();

    static CircularEdgeSPtr create(CircularVertexSPtr src, CircularVertexSPtr dst);

    CircularVertexSPtr getVertexSrc() const;
    void setVertexSrc(CircularVertexSPtr src);
    CircularVertexSPtr getVertexDst() const;
    void setVertexDst(CircularVertexSPtr dst);
    SphericalPolygonSPtr getPolygon() const;
    void setPolygon(SphericalPolygonSPtr polygon);
    std::list<CircularEdgeSPtr>::iterator getListIt() const;
    void setListIt(std::list<CircularEdgeSPtr>::iterator list_it);
    CircularEdgeDataSPtr getData() const;
    void setData(CircularEdgeDataSPtr data);
    bool hasData() const;

    Plane3SPtr getSupportingPlane() const;
    void setSupportingPlane(Plane3SPtr supporting_plane);
    bool initSupportingPlane();
    Plane3SPtr supportingPlane();

    CircularEdgeSPtr next() const;
    CircularEdgeSPtr prev() const;

    std::string toString() const;

protected:
    CircularEdge(CircularVertexSPtr src, CircularVertexSPtr dst);

    CircularVertexSPtr vertex_src_;
    CircularVertexSPtr vertex_dst_;
    SphericalPolygonWPtr polygon_;
    std::list<CircularEdgeSPtr>::iterator list_it_;
    CircularEdgeDataSPtr data_;
    Plane3SPtr supporting_plane_;
};

} }

#endif /* DATA_3D_CIRCULAREDGE_H */

