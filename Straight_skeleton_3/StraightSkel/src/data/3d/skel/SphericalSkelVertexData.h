/**
 * @file   data/3d/skel/SphericalSkelVertexData.h
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#ifndef DATA_3D_SKEL_SPHERICALSKELVERTEXDATA_H
#define DATA_3D_SKEL_SPHERICALSKELVERTEXDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/CircularVertexData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SphericalSkelVertexData : public CircularVertexData {
public:
    virtual ~SphericalSkelVertexData();

    static SphericalSkelVertexDataSPtr create(CircularVertexSPtr vertex);

    CircularArcSPtr getArc() const;
    void setArc(CircularArcSPtr arc);
    CircularNodeSPtr getNode() const;
    void setNode(CircularNodeSPtr node);
    CircularVertexSPtr getOffsetVertex() const;
    void setOffsetVertex(CircularVertexSPtr offset_vertex);

    double getSpeed() const;
    void setSpeed(double speed);

    EdgeSPtr getEdgeOrigin() const;
    void setEdgeOrigin(EdgeSPtr edge_origin);

    bool isReflex() const;
    void setReflex(bool reflex);

    bool doInvert() const;
    void setInvert(bool invert);

protected:
    SphericalSkelVertexData();
    CircularArcWPtr arc_;
    CircularNodeWPtr node_;
    CircularVertexWPtr offset_vertex_;
    double speed_;
    EdgeWPtr edge_origin_;
    bool reflex_;
    bool invert_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALSKELVERTEXDATA_H */

