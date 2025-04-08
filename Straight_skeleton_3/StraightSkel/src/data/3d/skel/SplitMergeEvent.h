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
 * @file   data/3d/skel/SplitMergeEvent.h
 * @author Gernot Walzl
 * @date   2013-08-09
 */

#ifndef DATA_3D_SKEL_SPLITMERGEEVENT_H
#define DATA_3D_SKEL_SPLITMERGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SplitMergeEvent : public AbstractEvent {
public:
    virtual ~SplitMergeEvent();
    static SplitMergeEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    VertexSPtr getVertex1() const;
    void setVertex1(VertexSPtr vertex_1);
    VertexSPtr getVertex2() const;
    void setVertex2(VertexSPtr vertex_2);
    FacetSPtr getFacet1() const;
    void setFacet1(FacetSPtr facet_1);
    FacetSPtr getFacet2() const;
    void setFacet2(FacetSPtr facet_2);
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
    bool operator==(const SplitMergeEvent& other) const;
protected:
    SplitMergeEvent();
    NodeSPtr node_;
    VertexWPtr vertex_1_;
    VertexWPtr vertex_2_;
    FacetWPtr facet_1_;
    FacetWPtr facet_2_;
    VertexFacetNeighborhood neighborhood_1_;
    VertexFacetNeighborhood neighborhood_2_;
};

} } }

#endif /* DATA_3D_SKEL_SPLITMERGEEVENT_H */

