/**
 * @file   data/3d/skel/AbstractEvent.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_ABSTRACTEVENT_H
#define DATA_3D_SKEL_ABSTRACTEVENT_H

#include <list>
#include <string>
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class AbstractEvent {
public:
    virtual ~AbstractEvent();

    PolyhedronSPtr getPolyhedronResult() const;
    void setPolyhedronResult(PolyhedronSPtr polyhedron);
    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<AbstractEventSPtr>::iterator getListIt() const;
    void setListIt(std::list<AbstractEventSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    virtual void setHighlight(bool highlight);
    virtual double getOffset() const = 0;  // abstract

    static const int CONST_OFFSET_EVENT = 1;
    static const int SAVE_OFFSET_EVENT = 2;

    /** 1 edge vanish event */
    static const int EDGE_EVENT = 3;

    /** 2 edge vanish event */
    static const int EDGE_MERGE_EVENT = 4;

    /** 3 edge vanish event */
    static const int TRIANGLE_EVENT = 5;

    /** 4 edge vanish event */
    static const int DBL_EDGE_MERGE_EVENT = 6;

    /** 5 edge vanish event */
    static const int DBL_TRIANGLE_EVENT = 7;

    /** 6 edge vanish event */
    static const int TETRAHEDRON_EVENT = 8;

    /** vertex-vertex contact event I */
    static const int VERTEX_EVENT = 9;

    /** vertex-vertex contact event II */
    static const int FLIP_VERTEX_EVENT = 10;

    /** vertex-edge contact event */
    static const int SURFACE_EVENT = 11;

    /** vertex-vertex-edge contact event I */
    static const int POLYHEDRON_SPLIT_EVENT = 12;

    /** vertex-vertex-edge contact event II */
    static const int SPLIT_MERGE_EVENT = 13;

    /** edge-edge contact event */
    static const int EDGE_SPLIT_EVENT = 14;

    /** vertex-facet contact event */
    static const int PIERCE_EVENT = 15;

    virtual int getType() const;

    virtual std::string toString() const;

protected:
    AbstractEvent();

    PolyhedronSPtr polyhedron_result_;
    StraightSkeletonWPtr skel_;
    std::list<AbstractEventSPtr>::iterator list_it_;
    int type_;
    int id_;
};

} } }

#endif /* DATA_3D_ABSTRACTEVENT_H */
