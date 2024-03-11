/**
 * @file   data/3d/skel/SphericalOffset.cpp
 * @author Gernot Walzl
 * @date   2012-12-21
 */

#include "data/3d/skel/SphericalOffset.h"

namespace data { namespace _3d { namespace skel {

SphericalOffset::SphericalOffset(double offset) {
    this->offset_ = offset;
    this->inf_jump_ = false;
}

SphericalOffset::~SphericalOffset() {
    // intentionally does nothing
}

SphericalOffsetSPtr SphericalOffset::create(double offset) {
    SphericalOffsetSPtr result = SphericalOffsetSPtr(new SphericalOffset(offset));
    return result;
}

double SphericalOffset::getOffset() const {
    return offset_;
}

void SphericalOffset::setOffset(double offset) {
    this->offset_ = offset;
}

bool SphericalOffset::isInfJump() const {
    return inf_jump_;
}

void SphericalOffset::setInfJump(bool inf_jump) {
    this->inf_jump_ = inf_jump;
}

int SphericalOffset::compareTo(SphericalOffsetSPtr other) {
    int result = 0;
    if (!this->inf_jump_ && other->inf_jump_) {
        result = -1;
    } else if (this->inf_jump_ && !other->inf_jump_) {
        result = 1;
    } else {
        if (this->offset_ < other->inf_jump_) {
            result = -1;
        } else if (this->offset_ > other->offset_) {
            result = 1;
        }
    }
    return result;
}

} } }
