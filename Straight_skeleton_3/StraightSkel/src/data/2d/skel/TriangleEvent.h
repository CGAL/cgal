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
 * @file   data/2d/skel/TriangleEvent.h
 * @author Gernot Walzl
 * @date   2013-05-06
 */

#ifndef DATA_2D_SKEL_TRIANGLEEVENT_H
#define DATA_2D_SKEL_TRIANGLEEVENT_H

#include "data/2d/skel/EdgeEvent.h"

namespace data { namespace _2d { namespace skel {

class TriangleEvent : public EdgeEvent {
public:
    virtual ~TriangleEvent();
    static TriangleEventSPtr create();
    void getVertices(VertexSPtr out[3]) const;
    void getEdges(EdgeSPtr out[3]) const;
    void setHighlight(bool highlight);
protected:
    TriangleEvent();
};

} } }

#endif /* DATA_2D_SKEL_TRIANGLEEVENT_H */
