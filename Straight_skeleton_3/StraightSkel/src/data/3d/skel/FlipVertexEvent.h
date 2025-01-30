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
 * @file   data/3d/skel/FlipVertexEvent.h
 * @author Gernot Walzl
 * @date   2012-11-22
 */

#ifndef DATA_3D_SKEL_FLIPVERTEXEVENT_H
#define DATA_3D_SKEL_FLIPVERTEXEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class FlipVertexEvent : public AbstractEvent {
public:
    virtual ~FlipVertexEvent();
    static FlipVertexEventSPtr create();
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
protected:
    FlipVertexEvent();
    NodeSPtr node_;
    VertexWPtr vertex_1_;
    VertexWPtr vertex_2_;
    FacetWPtr facet_1_;
    FacetWPtr facet_2_;
};

} } }

#endif /* DATA_3D_SKEL_FLIPVERTEXEVENT_H */

