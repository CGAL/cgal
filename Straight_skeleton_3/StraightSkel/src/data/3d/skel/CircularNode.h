/**
 * @file   data/3d/skel/CircularNode.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_SKEL_CIRCULARNODE_H
#define DATA_3D_SKEL_CIRCULARNODE_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class CircularNode : public std::enable_shared_from_this<CircularNode> {
public:
    virtual ~CircularNode();

    static CircularNodeSPtr create(Point3SPtr point);

    Point3SPtr getPoint() const;
    void setPoint(Point3SPtr point);

    double getOffset() const;
    void setOffset(double offset);

    SphericalSkeletonSPtr getSkel() const;
    void setSkel(SphericalSkeletonSPtr skel);
    std::list<CircularNodeSPtr>::iterator getListIt() const;
    void setListIt(std::list<CircularNodeSPtr>::iterator list_it);

    void addArc(CircularArcSPtr arc);
    bool removeArc(CircularArcSPtr arc);
    bool containsArc(CircularArcSPtr arc) const;

    void clear();

    std::list<CircularArcWPtr>& arcs();

    unsigned int degree();

    double getX() const;
    double getY() const;
    double getZ() const;

    std::string toString() const;

protected:
    CircularNode(Point3SPtr point);
    Point3SPtr point_;
    double offset_;
    std::list<CircularArcWPtr> arcs_;  // every CircularNode has 3 CircularArcs
    SphericalSkeletonWPtr skel_;
    std::list<CircularNodeSPtr>::iterator list_it_;
};

} } }

#endif /* DATA_3D_SKEL_CIRCULARNODE_H */
