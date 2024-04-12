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
 * @file   data/3d/skel/SaveOffsetEvent.cpp
 * @author Gernot Walzl
 * @date   2023-12-28
 */

#include "data/3d/skel/SaveOffsetEvent.h"

namespace data { namespace _3d { namespace skel {

SaveOffsetEvent::SaveOffsetEvent() {
    this->type_ = AbstractEvent::SAVE_OFFSET_EVENT;
    this->offset_ = -1.0;
}

SaveOffsetEvent::SaveOffsetEvent(CGAL::FT offset) {
    this->type_ = AbstractEvent::SAVE_OFFSET_EVENT;
    this->offset_ = offset;
}

SaveOffsetEvent::~SaveOffsetEvent() {
    // intentionally does nothing
}

SaveOffsetEventSPtr SaveOffsetEvent::create() {
    SaveOffsetEventSPtr result = SaveOffsetEventSPtr(new SaveOffsetEvent());
    return result;
}

SaveOffsetEventSPtr SaveOffsetEvent::create(CGAL::FT offset) {
    SaveOffsetEventSPtr result = SaveOffsetEventSPtr(
            new SaveOffsetEvent(offset));
    return result;
}

CGAL::FT SaveOffsetEvent::getOffset() const {
    return this->offset_;
}

void SaveOffsetEvent::setOffset(CGAL::FT offset) {
    this->offset_ = offset;
}

} } }
