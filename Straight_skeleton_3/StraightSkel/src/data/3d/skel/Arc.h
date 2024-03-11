/**
 * @file   data/3d/skel/Arc.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_ARC_H
#define DATA_3D_SKEL_ARC_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class Arc {
public:
    virtual ~Arc();
    static ArcSPtr create(NodeSPtr node_src, Vector3SPtr direction);
    static ArcSPtr create(NodeSPtr node_src, NodeSPtr node_dst);

    NodeSPtr getNodeSrc() const;
    void setNodeSrc(NodeSPtr node_src);
    std::list<ArcWPtr>::iterator getNodeSrcListIt() const;
    void setNodeSrcListIt(std::list<ArcWPtr>::iterator node_src_list_it);
    NodeSPtr getNodeDst() const;
    void setNodeDst(NodeSPtr node_dst);
    std::list<ArcWPtr>::iterator getNodeDstListIt() const;
    void setNodeDstListIt(std::list<ArcWPtr>::iterator node_dst_list_it);
    Vector3SPtr getDirection() const;
    void setDirection(Vector3SPtr direction);
    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<ArcSPtr>::iterator getListIt() const;
    void setListIt(std::list<ArcSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    void addSheet(SheetSPtr sheet);
    bool removeSheet(SheetSPtr sheet);

    std::list<SheetWPtr>& sheets();

    Line3SPtr line() const;

    bool hasNodeDst() const;

    std::string toString() const;

protected:
    Arc(NodeSPtr node_src, Vector3SPtr direction);
    Arc(NodeSPtr node_src, NodeSPtr node_dst);
    NodeSPtr node_src_;
    std::list<ArcWPtr>::iterator node_src_list_it_;
    NodeSPtr node_dst_;
    std::list<ArcWPtr>::iterator node_dst_list_it_;
    Vector3SPtr direction_;
    std::list<SheetWPtr> sheets_; // every arc has 3 sheets
    StraightSkeletonWPtr skel_;
    std::list<ArcSPtr>::iterator list_it_;
    int id_;
};

} } }

#endif /* DATA_3D_SKEL_ARC_H */

