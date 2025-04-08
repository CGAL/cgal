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
 * @file   data/3d/skel/ConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_CONSTOFFSETEVENT_H
#define DATA_3D_SKEL_CONSTOFFSETEVENT_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class ConstOffsetEvent : public AbstractEvent {
public:
    virtual ~ConstOffsetEvent();
    static ConstOffsetEventSPtr create();
    static ConstOffsetEventSPtr create(CGAL::FT offset);
    CGAL::FT getOffset() const override;
    void setOffset(CGAL::FT offset);
    bool operator==(const ConstOffsetEvent& other) const;
protected:
    ConstOffsetEvent();
    ConstOffsetEvent(CGAL::FT offset);
    CGAL::FT offset_;
};

} } }

#endif /* DATA_3D_SKEL_CONSTOFFSETEVENT_H */
