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
 * @file   data/3d/skel/EdgeSplitEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_EDGESPLITEVENT_H
#define DATA_3D_SKEL_EDGESPLITEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class EdgeSplitEvent : public AbstractEvent {
public:
    virtual ~EdgeSplitEvent();
    static EdgeSplitEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const;
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    int getEdgeOrientation() const;
    void setEdgeOrientation(int edge_orientation);
    void setHighlight(bool highlight);
    std::string toString() const override;
    bool isValid() const override;
protected:
    EdgeSplitEvent();
    NodeSPtr node_;
    EdgeWPtr edge1_;
    EdgeWPtr edge2_;
    int edge_orientation_;
};

} } }

#endif /* DATA_3D_SKEL_EDGESPLITEVENT_H */

