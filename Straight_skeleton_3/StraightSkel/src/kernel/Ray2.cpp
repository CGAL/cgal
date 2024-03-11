/**
 * @file   kernel/Ray2.cpp
 * @author Gernot Walzl
 * @date   2012-02-13
 */

#include "kernel/Ray2.h"

namespace kernel {

Ray2::Ray2() {
    this->point_ = new Point2();
    this->direction_ = new Vector2();
}

Ray2::Ray2(const Point2& point, const Vector2& direction) {
    this->point_ = new Point2(point);
    this->direction_ = new Vector2(direction);
}

Ray2::Ray2(const Ray2& orig) {
    this->point_ = new Point2(*(orig.point_));
    this->direction_ = new Vector2(*(orig.direction_));
}

Ray2::~Ray2() {
    delete this->point_;
    delete this->direction_;
}

}
