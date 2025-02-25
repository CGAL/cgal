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
 * @file   data/3d/skel/TriangleEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_TRIANGLEEVENT_H
#define DATA_3D_TRIANGLEEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class TriangleEvent : public AbstractEvent {
public:
    virtual ~TriangleEvent();
    static TriangleEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    FacetSPtr getFacet() const;
    void setFacet(FacetSPtr facet);
    EdgeSPtr getEdgeBegin() const;
    void setEdgeBegin(EdgeSPtr edge_begin);
    void getVertices(VertexSPtr out[3]) const;
    void getEdges(EdgeSPtr out[3]) const;
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
protected:
    TriangleEvent();
    NodeSPtr node_;
    FacetWPtr facet_; // @todo shouldn't be needed, edge_begin_->getFacetL is enough
    EdgeWPtr edge_begin_;
    EdgeFacetNeighborhood neighborhood_; // this covers the four faces involved in the triangle event
};

} } }

#endif /* DATA_3D_SKEL_TRIANGLEEVENT_H */

