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
 * @file   data/3d/skel/SphericalConstOffsetEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#include "data/3d/skel/SphericalConstOffsetEvent.h"

namespace data { namespace _3d { namespace skel {

SphericalConstOffsetEvent::SphericalConstOffsetEvent(CGAL::FT offset) {
    this->type_ = SphericalAbstractEvent::CONST_OFFSET_EVENT;
    this->offset_ = offset;
}

SphericalConstOffsetEvent::~SphericalConstOffsetEvent() {
    // intentionally does nothing
}

SphericalConstOffsetEventSPtr SphericalConstOffsetEvent::create(CGAL::FT offset) {
    SphericalConstOffsetEventSPtr result = SphericalConstOffsetEventSPtr(
            new SphericalConstOffsetEvent(offset));
    return result;
}

CGAL::FT SphericalConstOffsetEvent::getOffset() const {
    return this->offset_;
}

void SphericalConstOffsetEvent::setOffset(CGAL::FT offset) {
    this->offset_ = offset;
}

} } }
