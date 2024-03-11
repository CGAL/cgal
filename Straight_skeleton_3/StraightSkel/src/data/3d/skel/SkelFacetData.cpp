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
    DEBUG_WPTR(offset_facet_);
    if (this->offset_facet_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->offset_facet_);
}

void SkelFacetData::setOffsetFacet(FacetSPtr offset_facet) {
    this->offset_facet_ = offset_facet;
}

FacetSPtr SkelFacetData::getFacetOrigin() const {
    DEBUG_WPTR(facet_origin_);
    if (this->facet_origin_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_origin_);
}

void SkelFacetData::setFacetOrigin(FacetSPtr facet_origin) {
    this->facet_origin_ = facet_origin;
}


double SkelFacetData::getSpeed() const {
    return speed_;
}

void SkelFacetData::setSpeed(double speed) {
    speed_ = speed;
}

} } }
