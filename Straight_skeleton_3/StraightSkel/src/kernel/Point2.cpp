/**
 * @file   kernel/Point2.cpp
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#include "kernel/Point2.h"

#include <stdexcept>

namespace kernel {

Point2::Point2() {
    this->x_ = 0.0;
    this->y_ = 0.0;
}

Point2::Point2(double x, double y) {
    this->x_ = x;
    this->y_ = y;
}

Point2::Point2(const Point2& orig) {
    this->x_ = orig.x_;
    this->y_ = orig.y_;
}

Point2::~Point2() {
    // intentionally does nothing
}

double Point2::getX(void) const {
    return this->x_;
}

double Point2::getY(void) const {
    return this->y_;
}

double Point2::operator[](unsigned int i) const {
    if (i == 0) {
        return this->x_;
    } else if (i == 1) {
        return this->y_;
    } else {
        throw std::out_of_range("Index out of bounds.");
    }
}

Vector2 Point2::operator-(const Point2& q) const {
    double v_x = this->x_ - q.x_;
    double v_y = this->y_ - q.y_;
    return Vector2(v_x, v_y);
}

Point2 Point2::operator+(const Vector2& v) const {
    double x = this->x_ + v[0];
    double y = this->y_ + v[1];
    return Point2(x, y);
}

Point2 Point2::operator-(const Vector2& v) const {
    double x = this->x_ - v[0];
    double y = this->y_ - v[1];
    return Point2(x, y);
}

bool Point2::operator==(const Point2& q) const {
    return (this->x_ == q.getX() && this->y_ == q.getY());
}

}
