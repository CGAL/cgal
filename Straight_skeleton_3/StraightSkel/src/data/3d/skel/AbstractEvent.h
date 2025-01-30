// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

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
    virtual CGAL::FT getOffset() const = 0; // abstract

    static const int SAVE_OFFSET_EVENT = 0;

    static const int CONST_OFFSET_EVENT = 1;

    /** 1 edge vanish event */
    static const int EDGE_EVENT = 2;

    /** 2 edge vanish event */
    static const int EDGE_MERGE_EVENT = 3;

    /** 3 edge vanish event */
    static const int TRIANGLE_EVENT = 4;

    /** 4 edge vanish event */
    static const int DBL_EDGE_MERGE_EVENT = 5;

    /** 5 edge vanish event */
    static const int DBL_TRIANGLE_EVENT = 6;

    /** 6 edge vanish event */
    static const int TETRAHEDRON_EVENT = 7;

    /** vertex-vertex contact event I */
    static const int VERTEX_EVENT = 8;

    /** vertex-vertex contact event II */
    static const int FLIP_VERTEX_EVENT = 9;

    /** vertex-edge contact event */
    static const int SURFACE_EVENT = 10;

    /** vertex-vertex-edge contact event I */
    static const int POLYHEDRON_SPLIT_EVENT = 11;

    /** vertex-vertex-edge contact event II */
    static const int SPLIT_MERGE_EVENT = 12;

    /** edge-edge contact event */
    static const int EDGE_SPLIT_EVENT = 13;

    /** vertex-facet contact event */
    static const int PIERCE_EVENT = 14;

    virtual int getType() const;

    virtual std::string toString() const;

    virtual bool isValid() const;

protected:
    AbstractEvent();

    PolyhedronSPtr polyhedron_result_;
    StraightSkeletonWPtr skel_;
    std::list<AbstractEventSPtr>::iterator list_it_;
    int type_;
    int id_;

private:
    static int next_id_;
};

class AbstractEventSPtrCompare
  : public CGAL::cpp98::binary_function<bool, AbstractEventSPtr, AbstractEventSPtr>
{
public:
    bool operator()(const AbstractEventSPtr& eventA,
                    const AbstractEventSPtr& eventB) const
    {
        // save and const events are handled outside of the queue, and in generic
        // position, all events have their own offsets
        CGAL_assertion(eventA->getOffset() != eventB->getOffset());
        return (eventA->getOffset() < eventB->getOffset());
    }
};

} } }

#endif /* DATA_3D_ABSTRACTEVENT_H */
