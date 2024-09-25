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
 * @file   data/3d/skel/DblEdgeMergeEvent.h
 * @author Gernot Walzl
 * @date   2012-10-30
 */

#ifndef DATA_3D_SKEL_DBLEDGEMERGEEVENT_H
#define DATA_3D_SKEL_DBLEDGEMERGEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class DblEdgeMergeEvent : public AbstractEvent {
public:
    virtual ~DblEdgeMergeEvent();
    static DblEdgeMergeEventSPtr create(PolyhedronSPtr polyhedron);
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const;
    FacetSPtr getFacet1() const;
    void setFacet1(FacetSPtr facet_1);
    EdgeSPtr getEdge11() const;
    void setEdge11(EdgeSPtr edge_11);
    EdgeSPtr getEdge12() const;
    void setEdge12(EdgeSPtr edge_12);
    FacetSPtr getFacet2() const;
    void setFacet2(FacetSPtr facet_2);
    EdgeSPtr getEdge21() const;
    void setEdge21(EdgeSPtr edge_21);
    EdgeSPtr getEdge22() const;
    void setEdge22(EdgeSPtr edge_22);
    void getVertices(VertexSPtr out[4]) const;
    void getEdges(EdgeSPtr out[4]) const;
    void setHighlight(bool highlight);
    std::string toString() const override;
protected:
    DblEdgeMergeEvent(PolyhedronSPtr polyhedron);
    NodeSPtr node_;
    FacetSPtr facet_1_;
    EdgeSPtr edge_11_;
    EdgeSPtr edge_12_;
    FacetSPtr facet_2_;
    EdgeSPtr edge_21_;
    EdgeSPtr edge_22_;
};

} } }

#endif /* DATA_3D_SKEL_DBLEDGEMERGEEVENT_H */

