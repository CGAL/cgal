// Copyright (c) 2024 GeometryFactory (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé

#ifndef DATA_3D_SKEL_VANISHEVENT_H
#define DATA_3D_SKEL_VANISHEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

#include <string>

namespace data { namespace _3d { namespace skel {

// @todo also add all the other boiler code as for other events (DAO, etc.)
class VanishEvent : public AbstractEvent {
public:
    virtual ~VanishEvent();
    static VanishEventSPtr create();
    NodeSPtr getNode() const;
    void setNode(NodeSPtr node);
    CGAL::FT getOffset() const override;
    EdgeSPtr getEdge() const;
    void setEdge(EdgeSPtr edge);
    void setHighlight(bool highlight) override;
    std::string toString() const override;
    bool isValid() const override;
    bool isObsolete() const override;
    bool operator==(const VanishEvent& other) const;
protected:
    VanishEvent();

    NodeSPtr node_;
    EdgeWPtr edge_;
    EdgeFacetNeighborhood neighborhood_;
};

} } }

#endif /* DATA_3D_SKEL_VANISHEVENT_H */

