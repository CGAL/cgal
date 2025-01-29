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
 * @file   data/3d/skel/ConstOffsetEvent.cpp
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#include "data/3d/skel/ConstOffsetEvent.h"

namespace data { namespace _3d { namespace skel {

ConstOffsetEvent::ConstOffsetEvent() {
    this->type_ = AbstractEvent::CONST_OFFSET_EVENT;
    this->offset_ = -1.0;
}

ConstOffsetEvent::ConstOffsetEvent(CGAL::FT offset) {
    this->type_ = AbstractEvent::CONST_OFFSET_EVENT;
    this->offset_ = offset;
}

ConstOffsetEvent::~ConstOffsetEvent() {
    // intentionally does nothing
}

ConstOffsetEventSPtr ConstOffsetEvent::create() {
    ConstOffsetEventSPtr result = ConstOffsetEventSPtr(new ConstOffsetEvent());
    return result;
}

ConstOffsetEventSPtr ConstOffsetEvent::create(CGAL::FT offset) {
    ConstOffsetEventSPtr result = ConstOffsetEventSPtr(
            new ConstOffsetEvent(offset));
    return result;
}

CGAL::FT ConstOffsetEvent::getOffset() const {
    return this->offset_;
}

void ConstOffsetEvent::setOffset(CGAL::FT offset) {
    this->offset_ = offset;
}

} } }
