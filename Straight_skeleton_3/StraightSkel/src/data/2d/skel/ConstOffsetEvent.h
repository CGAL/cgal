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
 * @file   data/2d/skel/ConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-04-06
 */

#ifndef DATA_2D_SKEL_CONSTOFFSETEVENT_H
#define DATA_2D_SKEL_CONSTOFFSETEVENT_H

#include "debug.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/AbstractEvent.h"

namespace data { namespace _2d { namespace skel {

class ConstOffsetEvent : public AbstractEvent {
public:
    virtual ~ConstOffsetEvent();

    static ConstOffsetEventSPtr create();
    static ConstOffsetEventSPtr create(CGAL::FT offset);

    CGAL::FT getOffset() const;
    void setOffset(CGAL::FT offset);
protected:
    ConstOffsetEvent();
    ConstOffsetEvent(CGAL::FT offset);
    CGAL::FT offset_;
};

} } }

#endif /* DATA_2D_SKEL_CONSTOFFSETEVENT_H */

