/**
 * @file   kernel/Segment3.cpp
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#include "kernel/Segment3.h"

namespace kernel {

Segment3::Segment3(const Point3& p, const Point3& q) {
    this->p_ = &p;
    this->q_ = &q;
}

Segment3::Segment3(const Segment3& orig) {
    this->p_ = orig.p_;
    this->q_ = orig.q_;
}

Segment3::~Segment3() {
    // intentionally does nothing
}

const Point3& Segment3::getP() const {
    return *(this->p_);
}

const Point3& Segment3::getQ() const {
    return *(this->q_);
}

void Segment3::setP(const Point3& p) {
    this->p_ = &p;
}

void Segment3::setQ(const Point3& q) {
    this->q_ = &q;
}

Line3 Segment3::line() const {
    return Line3(*(this->p_), *(this->q_));
}

}
