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
 * @file   data/3d/CircularEdgeData.h
 * @author Gernot Walzl
 * @date   2012-11-30
 */

#ifndef DATA_3D_CIRCULAREDGEDATA_H
#define DATA_3D_CIRCULAREDGEDATA_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

/*!
 * This class is intended to be subclassed.
 */
class CircularEdgeData {
public:
    virtual ~CircularEdgeData();

    CircularEdgeSPtr getEdge() const;
    void setEdge(CircularEdgeSPtr edge);

    bool isHighlight() const;
    void setHighlight(bool highlight);

protected:
    CircularEdgeData();
    CircularEdgeWPtr edge_;
    bool highlight_;
};

} }

#endif /* DATA_3D_CIRCULAREDGEDATA_H */

