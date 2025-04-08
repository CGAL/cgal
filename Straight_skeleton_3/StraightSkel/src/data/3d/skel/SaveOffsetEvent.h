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
 * @file   data/3d/skel/SaveOffsetEvent.h
 * @author Gernot Walzl
 * @date   2013-12-28
 */

#ifndef DATA_3D_SKEL_SAVEOFFSETEVENT_H
#define DATA_3D_SKEL_SAVEOFFSETEVENT_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SaveOffsetEvent : public AbstractEvent {
public:
    virtual ~SaveOffsetEvent();
    static SaveOffsetEventSPtr create();
    static SaveOffsetEventSPtr create(CGAL::FT offset);
    CGAL::FT getOffset() const override;
    void setOffset(CGAL::FT offset);
    bool operator==(const SaveOffsetEvent& other) const;
protected:
    SaveOffsetEvent();
    SaveOffsetEvent(CGAL::FT offset);
    CGAL::FT offset_;
};

} } }

#endif /* DATA_3D_SKEL_SAVEOFFSETEVENT_H */
