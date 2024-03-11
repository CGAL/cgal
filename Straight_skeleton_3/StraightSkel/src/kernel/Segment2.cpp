/**
 * @file   kernel/Segment2.cpp
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#include "kernel/Segment2.h"

namespace kernel {

Segment2::Segment2(const Point2& p, const Point2& q) {
    this->p_ = &p;
    this->q_ = &q;
}

Segment2::Segment2(const Segment2& orig) {
    this->p_ = orig.p_;
    this->q_ = orig.q_;
}

Segment2::~Segment2() {
    // intentionally does nothing
}

const Point2& Segment2::getP() const {
    return *(this->p_);
}

const Point2& Segment2::getQ() const {
    return *(this->q_);
}

void Segment2::setP(const Point2& p) {
    this->p_ = &p;
}

void Segment2::setQ(const Point2& q) {
    this->q_ = &q;
}

Line2 Segment2::line() const {
    return Line2(*(this->p_), *(this->q_));
}

}
