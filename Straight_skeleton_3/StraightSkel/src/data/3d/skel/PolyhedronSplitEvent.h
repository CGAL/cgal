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
 * @file   data/3d/skel/PolyhedronSplitEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_POLYHEDRONSPLITEVENT_H
#define DATA_3D_SKEL_POLYHEDRONSPLITEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class PolyhedronSplitEvent : public AbstractEvent {
public:
    virtual ~PolyhedronSplitEvent();
    static PolyhedronSplitEventSPtr create(PolyhedronSPtr polyhedron);
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const;
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    void setHighlight(bool highlight);
    std::string toString() const override;
protected:
    PolyhedronSplitEvent(PolyhedronSPtr polyhedron);
    NodeSPtr node_;
    EdgeSPtr edge1_;
    EdgeSPtr edge2_;
};

} } }

#endif /* DATA_3D_SKEL_POLYHEDRONSPLITEVENT_H */

