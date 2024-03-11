/**
 * @file   data/2d/skel/StraightSkeleton.h
 * @author Gernot Walzl
 * @date   2011-11-21
 */

#ifndef DATA_2D_SKEL_STRAIGHTSKELETON_H
#define DATA_2D_SKEL_STRAIGHTSKELETON_H

#include "typedefs_thread.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace skel {

class StraightSkeleton : public std::enable_shared_from_this<StraightSkeleton> {
public:
    virtual ~StraightSkeleton();

    static StraightSkeletonSPtr create();

    /**
     * also adds the node
     */
    void addEvent(AbstractEventSPtr event);
    bool removeEvent(AbstractEventSPtr event);
    void addNode(NodeSPtr node);
    bool removeNode(NodeSPtr node);
    void addArc(ArcSPtr arc);
    bool removeArc(ArcSPtr arc);

    SharedMutex& mutex();
    std::list<AbstractEventSPtr>& events();
    std::list<NodeSPtr>& nodes();
    std::list<ArcSPtr>& arcs();

    PolygonSPtr getPolygon() const;
    void setPolygon(PolygonSPtr polygon);

    int getID() const;
    void setID(int id);

    void resetAllIDs();

    bool isConsistent() const;

    int countEvents(int type) const;

    std::string getDescription() const;
    void setDescription(const std::string& desc);
    void appendDescription(const std::string& desc);

    std::string toString() const;

protected:
    StraightSkeleton();
    PolygonSPtr polygon_;
    mutable SharedMutex mutex_;
    std::list<AbstractEventSPtr> events_;
    std::list<NodeSPtr> nodes_;
    std::list<ArcSPtr> arcs_;
    int id_;
    std::string description_;
};

} } }

#endif /* DATA_2D_SKEL_STRAIGHTSKELETON_H */
