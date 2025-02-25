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
 * @file   data/3d/skel/EdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2012-09-14
 */

#ifndef DATA_3D_SKEL_EDGEMERGEEVENT_H
#define DATA_3D_SKEL_EDGEMERGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class EdgeMergeEvent : public AbstractEvent {
public:
    virtual ~EdgeMergeEvent();
    static EdgeMergeEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    FacetSPtr getFacet() const;
    void setFacet(FacetSPtr facet);
    EdgeSPtr getEdge1() const;
    void setEdge1(EdgeSPtr edge1);
    EdgeSPtr getEdge2() const;
    void setEdge2(EdgeSPtr edge2);
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
protected:
    EdgeMergeEvent();
    NodeSPtr node_;
    FacetWPtr facet_;
    EdgeWPtr edge1_;
    EdgeWPtr edge2_;
};

} } }

#endif /* DATA_3D_SKEL_EDGEMERGEEVENT_H */

