/**
 * @file   data/3d/skel/SphericalConstOffsetEvent.cpp
 * @author Gernot Walzl
 * @date   2012-11-29
 */

#include "data/3d/skel/SphericalConstOffsetEvent.h"

namespace data { namespace _3d { namespace skel {

SphericalConstOffsetEvent::SphericalConstOffsetEvent(double offset) {
    this->type_ = SphericalAbstractEvent::CONST_OFFSET_EVENT;
    this->offset_ = offset;
}

SphericalConstOffsetEvent::~SphericalConstOffsetEvent() {
    // intentionally does nothing
}

SphericalConstOffsetEventSPtr SphericalConstOffsetEvent::create(double offset) {
    SphericalConstOffsetEventSPtr result = SphericalConstOffsetEventSPtr(
            new SphericalConstOffsetEvent(offset));
    return result;
}

double SphericalConstOffsetEvent::getOffset() const {
    return this->offset_;
}

void SphericalConstOffsetEvent::setOffset(double offset) {
    this->offset_ = offset;
}

} } }
