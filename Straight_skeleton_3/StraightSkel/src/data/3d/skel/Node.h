/**
 * @file   data/3d/skel/Node.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_NODE_H
#define DATA_3D_SKEL_NODE_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class Node : public std::enable_shared_from_this<Node> {
public:
    virtual ~Node();
    static NodeSPtr create(Point3SPtr point);

    Point3SPtr getPoint() const;
    void setPoint(Point3SPtr point);
    double getOffset() const;
    void setOffset(double offset);

    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<NodeSPtr>::iterator getListIt() const;
    void setListIt(std::list<NodeSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    void addArc(ArcSPtr arc);
    bool removeArc(ArcSPtr arc);

    void addSheet(SheetSPtr sheet);
    bool removeSheet(SheetSPtr sheet);

    bool containsArc(ArcSPtr arc) const;
    bool containsSheet(SheetSPtr sheet) const;

    void clear();

    std::list<ArcWPtr>& arcs();
    std::list<SheetWPtr>& sheets();

    unsigned int degree() const;

    double getX() const;
    double getY() const;
    double getZ() const;

    std::string toString() const;

protected:
    Node(Point3SPtr point);
    Point3SPtr point_;
    double offset_;
    std::list<ArcWPtr> arcs_;
    std::list<SheetWPtr> sheets_;
    StraightSkeletonWPtr skel_;
    std::list<NodeSPtr>::iterator list_it_;
    int id_;
};

} } }

#endif /* DATA_3D_SKEL_NODE_H */
