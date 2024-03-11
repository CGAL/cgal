/**
 * @file   data/3d/skel/SphericalInversionEvent.cpp
 * @author Gernot Walzl
 * @date   2013-03-13
 */

#include "data/3d/skel/SphericalInversionEvent.h"

namespace data { namespace _3d { namespace skel {

SphericalInversionEvent::SphericalInversionEvent() {
    this->type_ = SphericalAbstractEvent::INVERSION_EVENT;
    this->offset_ = 0.0;
}

SphericalInversionEvent::~SphericalInversionEvent() {
    polygon_.reset();
}

SphericalInversionEventSPtr SphericalInversionEvent::create() {
    SphericalInversionEventSPtr result = SphericalInversionEventSPtr(
            new SphericalInversionEvent());
    return result;
}

double SphericalInversionEvent::getOffset() const {
    return this->offset_;
}

void SphericalInversionEvent::setOffset(double offset) {
    this->offset_ = offset;
}

SphericalPolygonSPtr SphericalInversionEvent::getPolygon() const {
    return this->polygon_;
}

void SphericalInversionEvent::setPolygon(SphericalPolygonSPtr polygon) {
    this->polygon_ = polygon;
}

} } }
