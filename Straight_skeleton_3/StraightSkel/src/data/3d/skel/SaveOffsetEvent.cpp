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

SaveOffsetEvent::SaveOffsetEvent(double offset) {
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

SaveOffsetEventSPtr SaveOffsetEvent::create(double offset) {
    SaveOffsetEventSPtr result = SaveOffsetEventSPtr(
            new SaveOffsetEvent(offset));
    return result;
}

double SaveOffsetEvent::getOffset() const {
    return this->offset_;
}

void SaveOffsetEvent::setOffset(double offset) {
    this->offset_ = offset;
}

} } }
