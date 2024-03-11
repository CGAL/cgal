/**
 * @file   data/3d/FacetData.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/FacetData.h"

#include "data/3d/Facet.h"

namespace data { namespace _3d {

FacetData::FacetData() {
    highlight_ = false;
}

FacetData::~FacetData() {
    // intentionally does nothing
}

FacetDataSPtr FacetData::create(FacetSPtr facet) {
    FacetDataSPtr result = FacetDataSPtr(new FacetData());
    result->setFacet(facet);
    facet->setData(result);
    return result;
}

FacetSPtr FacetData::getFacet() const {
    if (this->facet_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_);
}

void FacetData::setFacet(FacetSPtr facet) {
    this->facet_ = facet;
}

bool FacetData::isHighlight() const {
    return highlight_;
}

void FacetData::setHighlight(bool highlight) {
    highlight_ = highlight;
}

} }
