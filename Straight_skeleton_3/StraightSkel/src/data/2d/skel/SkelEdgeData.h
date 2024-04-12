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
 * @file   data/2d/skel/SkelEdgeData.h
 * @author Gernot Walzl
 * @date   2012-05-12
 */

#ifndef DATA_2D_SKEL_SKELEDGEDATA_H
#define DATA_2D_SKEL_SKELEDGEDATA_H

#include "data/2d/ptrs.h"
#include "data/2d/EdgeData.h"
#include "data/2d/skel/ptrs.h"

namespace data { namespace _2d { namespace skel {

class SkelEdgeData : public EdgeData {
public:
    virtual ~SkelEdgeData();

    static SkelEdgeDataSPtr create(EdgeSPtr edge);

    EdgeSPtr getOffsetEdge() const;
    void setOffsetEdge(EdgeSPtr offset_edge);

    EdgeSPtr getEdgeOrigin() const;
    void setEdgeOrigin(EdgeSPtr edge_origin);

    CGAL::FT getSpeed() const;
    void setSpeed(CGAL::FT speed);

protected:
    SkelEdgeData();
    EdgeWPtr offset_edge_;
    EdgeWPtr edge_origin_;
    CGAL::FT speed_;
};

} } }

#endif /* DATA_2D_SKEL_SKELEDGEDATA_H */

