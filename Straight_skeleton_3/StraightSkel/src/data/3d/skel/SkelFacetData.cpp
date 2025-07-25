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
 * @file   data/3d/skel/SkelFacetData.cpp
 * @author Gernot Walzl
 * @date   2012-05-07
 */

#include "data/3d/skel/SkelFacetData.h"

#include "debug.h"
#include "data/3d/Facet.h"

namespace data { namespace _3d { namespace skel {

SkelFacetData::SkelFacetData() {
    speed_ = 1.0;
}

SkelFacetData::~SkelFacetData() {
    // intentionally does nothing
}

SkelFacetDataSPtr SkelFacetData::create(FacetSPtr facet) {
    SkelFacetDataSPtr result = SkelFacetDataSPtr(new SkelFacetData());
    result->setFacet(facet);
    result->setFacetOrigin(facet);
    facet->setData(result);
    return result;
}

FacetSPtr SkelFacetData::getOffsetFacet() const {
    CGAL_SS3_DEBUG_WPTR(offset_facet_);
    if (this->offset_facet_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->offset_facet_);
}

void SkelFacetData::setOffsetFacet(FacetSPtr offset_facet) {
    this->offset_facet_ = offset_facet;
}

FacetSPtr SkelFacetData::getFacetOrigin() const {
    CGAL_SS3_DEBUG_WPTR(facet_origin_);
    if (this->facet_origin_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_origin_);
}

void SkelFacetData::setFacetOrigin(FacetSPtr facet_origin) {
    this->facet_origin_ = facet_origin;
}


CGAL::FT SkelFacetData::getSpeed() const {
    CGAL_assertion(speed_ != 0);
    return speed_;
}

void SkelFacetData::setSpeed(CGAL::FT speed) {
    speed_ = speed;
}

int SkelFacetData::getStepID() const {
  return this->step_id_;
}

void SkelFacetData::setStepID(int step_id) {
  this->step_id_ = step_id;
}

} } }
