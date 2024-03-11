/**
 * @file   data/3d/skel/CircularArc.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_SKEL_CIRCULARARC_H
#define DATA_3D_SKEL_CIRCULARARC_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class CircularArc {
public:
    virtual ~CircularArc();

    static CircularArcSPtr create(CircularNodeSPtr node_src, Vector3SPtr direction);
    static CircularArcSPtr create(CircularNodeSPtr node_src, CircularNodeSPtr node_dst);

    CircularNodeSPtr getNodeSrc() const;
    void setNodeSrc(CircularNodeSPtr node_src);
    std::list<CircularArcWPtr>::iterator getNodeSrcListIt() const;
    void setNodeSrcListIt(std::list<CircularArcWPtr>::iterator node_src_list_it);
    CircularNodeSPtr getNodeDst() const;
    void setNodeDst(CircularNodeSPtr node_dst);
    std::list<CircularArcWPtr>::iterator getNodeDstListIt() const;
    void setNodeDstListIt(std::list<CircularArcWPtr>::iterator node_dst_list_it);
    Vector3SPtr getDirection() const;
    void setDirection(Vector3SPtr direction);
    CircularEdgeSPtr getEdgeLeft() const;
    void setEdgeLeft(CircularEdgeSPtr edge_left);
    CircularEdgeSPtr getEdgeRight() const;
    void setEdgeRight(CircularEdgeSPtr edge_right);
    Plane3SPtr getSupportingPlane() const;
    void setSupportingPlane(Plane3SPtr supporting_plane);
    SphericalSkeletonSPtr getSkel() const;
    void setSkel(SphericalSkeletonSPtr skel);
    std::list<CircularArcSPtr>::iterator getListIt() const;
    void setListIt(std::list<CircularArcSPtr>::iterator list_it);

    bool hasNodeDst() const;

    std::string toString() const;

protected:
    CircularArc(CircularNodeSPtr node_src, Vector3SPtr direction);
    CircularArc(CircularNodeSPtr node_src, CircularNodeSPtr node_dst);
    CircularNodeSPtr node_src_;
    std::list<CircularArcWPtr>::iterator node_src_list_it_;
    CircularNodeSPtr node_dst_;
    std::list<CircularArcWPtr>::iterator node_dst_list_it_;
    Vector3SPtr direction_;
    CircularEdgeSPtr edge_left_;
    CircularEdgeSPtr edge_right_;
    Plane3SPtr supporting_plane_;
    SphericalSkeletonWPtr skel_;
    std::list<CircularArcSPtr>::iterator list_it_;
};

} } }

#endif /* DATA_3D_SKEL_CIRCULARARC_H */

