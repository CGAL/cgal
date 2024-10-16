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
 * @file   data/3d/skel/EdgeEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_EDGEEVENT_H
#define DATA_3D_SKEL_EDGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class EdgeEvent : public AbstractEvent {
public:
    virtual ~EdgeEvent();
    static EdgeEventSPtr create(PolyhedronSPtr polyhedron);
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const;
    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);
    void setHighlight(bool highlight);
    std::string toString() const override;
protected:
    EdgeEvent(PolyhedronSPtr polyhedron);

    NodeSPtr node_;
    EdgeSPtr edge_;
};

} } }

#endif /* DATA_3D_SKEL_EDGEEVENT_H */

