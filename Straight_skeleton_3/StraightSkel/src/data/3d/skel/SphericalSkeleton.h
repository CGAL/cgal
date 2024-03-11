/**
 * @file   data/3d/skel/SphericalSkeleton.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_SKEL_SPHERICALSKELETON_H
#define DATA_3D_SKEL_SPHERICALSKELETON_H

#include "typedefs_thread.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class SphericalSkeleton : public std::enable_shared_from_this<SphericalSkeleton> {
public:
    virtual ~SphericalSkeleton();

    static SphericalSkeletonSPtr create(Sphere3SPtr sphere);

    Sphere3SPtr getSphere() const;

    Plane3SPtr getRotAxesPlane() const;
    void setRotAxesPlane(Plane3SPtr rot_axes_plane);

    /**
     * also adds the node
     */
    void addEvent(SphericalAbstractEventSPtr event);
    bool removeEvent(SphericalAbstractEventSPtr event);
    void addNode(CircularNodeSPtr node);
    bool removeNode(CircularNodeSPtr node);
    void addArc(CircularArcSPtr arc);
    bool removeArc(CircularArcSPtr arc);

    SharedMutex& mutex();
    std::list<SphericalAbstractEventSPtr>& events();
    std::list<CircularNodeSPtr>& nodes();
    std::list<CircularArcSPtr>& arcs();

    bool isConsistent() const;

    int countEvents(int type) const;

    double getRadius() const;

    std::string toString() const;

protected:
    SphericalSkeleton(Sphere3SPtr sphere);
    Sphere3SPtr sphere_;
    Plane3SPtr rot_axes_plane_;
    mutable SharedMutex mutex_;
    std::list<SphericalAbstractEventSPtr> events_;
    std::list<CircularNodeSPtr> nodes_;
    std::list<CircularArcSPtr> arcs_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALSKELETON_H */
