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
 * @file   data/3d/skel/TetrahedronEvent.h
 * @author Gernot Walzl
 * @date   2012-04-23
 */

#ifndef DATA_3D_SKEL_TETRAHEDRONEVENT_H
#define DATA_3D_SKEL_TETRAHEDRONEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

class TetrahedronEvent : public AbstractEvent {
public:
    virtual ~TetrahedronEvent();
    static TetrahedronEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    EdgeSPtr getEdgeBegin() const;
    void setEdgeBegin(EdgeSPtr edge_begin);
    void getVertices(VertexSPtr out[4]) const;
    void getEdges(EdgeSPtr out[6]) const;
    void getFacets(FacetSPtr out[4]) const;
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
protected:
    TetrahedronEvent();
    NodeSPtr node_;
    EdgeWPtr edge_begin_;
    EdgeFacetNeighborhood neighborhood_;
};

} } }

#endif /* DATA_3D_SKEL_TETRAHEDRONEVENT_H */

