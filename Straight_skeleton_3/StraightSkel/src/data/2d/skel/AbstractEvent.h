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
 * @file   data/2d/skel/AbstractEvent.h
 * @author Gernot Walzl
 * @date   2012-02-02
 */

#ifndef DATA_2D_SKEL_ABSTRACTEVENT_H
#define DATA_2D_SKEL_ABSTRACTEVENT_H

#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _2d { namespace skel {

class AbstractEvent {
public:
    virtual ~AbstractEvent();
    PolygonSPtr getPolygonResult() const;
    void setPolygonResult(PolygonSPtr polyon);
    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<AbstractEventSPtr>::iterator getListIt() const;
    void setListIt(std::list<AbstractEventSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    virtual void setHighlight(bool highlight);
    virtual CGAL::FT getOffset() const = 0;  // abstract

    static const int CONST_OFFSET_EVENT = 1;
    static const int EDGE_EVENT = 2;
    static const int SPLIT_EVENT = 3;
    static const int TRIANGLE_EVENT = 4;

    virtual int getType() const;

    virtual std::string toString() const;

protected:
    AbstractEvent();
    PolygonSPtr polygon_result_;
    StraightSkeletonWPtr skel_;
    std::list<AbstractEventSPtr>::iterator list_it_;
    int type_;
    int id_;
};

} } }

#endif /* DATA_2D_SKEL_ABSTRACTEVENT_H */

