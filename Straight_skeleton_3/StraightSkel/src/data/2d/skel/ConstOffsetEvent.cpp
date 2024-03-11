/**
 * @file   data/2d/skel/ConstOffsetEvent.cpp
 * @author Gernot Walzl
 * @date   2012-04-06
 */

#include "data/2d/skel/ConstOffsetEvent.h"

namespace data { namespace _2d { namespace skel {

ConstOffsetEvent::ConstOffsetEvent() {
    this->type_ = AbstractEvent::CONST_OFFSET_EVENT;
    this->offset_ = 1.0;
}

ConstOffsetEvent::ConstOffsetEvent(double offset) {
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

ConstOffsetEventSPtr ConstOffsetEvent::create(double offset) {
    ConstOffsetEventSPtr result = ConstOffsetEventSPtr(
            new ConstOffsetEvent(offset));
    return result;
}

double ConstOffsetEvent::getOffset() const {
    return this->offset_;
}

void ConstOffsetEvent::setOffset(double offset) {
    this->offset_ = offset;
}

} } }
