/**
 * @file   data/3d/skel/StraightSkeleton.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_SKEL_STRAIGHTSKELETON_H
#define DATA_3D_SKEL_STRAIGHTSKELETON_H

#include "typedefs_thread.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

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
    void addSheet(SheetSPtr sheet);
    bool removeSheet(SheetSPtr sheet);

    SharedMutex& mutex();
    std::list<AbstractEventSPtr>& events();
    std::list<NodeSPtr>& nodes();
    std::list<ArcSPtr>& arcs();
    std::list<SheetSPtr>& sheets();

    PolyhedronSPtr getPolyhedron() const;
    void setPolyhedron(PolyhedronSPtr polyhedron);

    int getID() const;
    void setID(int id);

    void resetAllIDs();

    std::string getConfig() const;
    void setConfig(const std::string& config);
    void appendConfig(const std::string& config);

    bool isConsistent() const;

    int countEvents(int type) const;

    std::string getDescription() const;
    void setDescription(const std::string& desc);
    void appendDescription(const std::string& desc);

    std::string toString() const;

protected:
    StraightSkeleton();
    PolyhedronSPtr polyhedron_;
    mutable SharedMutex mutex_;
    std::list<AbstractEventSPtr> events_;
    std::list<NodeSPtr> nodes_;
    std::list<ArcSPtr> arcs_;
    std::list<SheetSPtr> sheets_;
    int id_;
    std::string config_;
    std::string description_;
};

} } }

#endif /* DATA_3D_SKEL_STRAIGHTSKELETON_H */
