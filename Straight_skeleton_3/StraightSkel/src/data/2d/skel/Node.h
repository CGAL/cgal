/**
 * @file   data/2d/skel/Node.h
 * @author Gernot Walzl
 * @date   2012-02-03
 */

#ifndef DATA_2D_SKEL_NODE_H
#define DATA_2D_SKEL_NODE_H

#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace skel {

class Node : public std::enable_shared_from_this<Node> {
public:
    virtual ~Node();
    static NodeSPtr create(Point2SPtr point);

    Point2SPtr getPoint() const;
    void setPoint(Point2SPtr point);
    double getHeight() const;
    void setHeight(double height);
    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<NodeSPtr>::iterator getListIt() const;
    void setListIt(std::list<NodeSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    void addArc(ArcSPtr arc);
    bool removeArc(ArcSPtr arc);

    std::list<ArcWPtr>& arcs();

    unsigned int degree() const;

    double getX() const;
    double getY() const;

    std::string toString() const;

protected:
    Node(Point2SPtr point);
    Point2SPtr point_;
    double height_;
    StraightSkeletonWPtr skel_;
    std::list<NodeSPtr>::iterator list_it_;
    std::list<ArcWPtr> arcs_;
    int id_;
};

} } }

#endif /* DATA_2D_SKEL_NODE_H */
