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
    static PolyhedronSplitEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
protected:
    PolyhedronSplitEvent();
    NodeSPtr node_;
    EdgeWPtr edge1_;
    EdgeWPtr edge2_;
    EdgeFacetNeighborhood neighborhood1_;
    EdgeFacetNeighborhood neighborhood2_;
};

} } }

#endif /* DATA_3D_SKEL_POLYHEDRONSPLITEVENT_H */

