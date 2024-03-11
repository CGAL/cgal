/**
 * @file   kernel/Line2.cpp
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#include "kernel/Line2.h"

namespace kernel {

Line2::Line2() {
    this->a_ = 0.0;
    this->b_ = 0.0;
    this->c_ = 0.0;
}

Line2::Line2(const Line2& orig) {
    this->a_ = orig.a_;
    this->b_ = orig.b_;
    this->c_ = orig.c_;
}

Line2::Line2(const Point2& p, const Point2& q) {
    a_ = p.getY() - q.getY();
    b_ = q.getX() - p.getX();
    c_ = - a_*p.getX() - b_*p.getY();
}

Line2::Line2(const Point2& p, const Vector2& direction) {
    a_ = -direction[1];
    b_ = direction[0];
    c_ = - a_*p.getX() - b_*p.getY();
}

Line2::Line2(double a, double b, double c) {
    a_ = a;
    b_ = b;
    c_ = c;
}

Line2::~Line2() {
    // intentionally does nothing
}

double Line2::getA() const {
    return this->a_;
}

double Line2::getB() const {
    return this->b_;
}

double Line2::getC() const {
    return this->c_;
}

Point2 Line2::point() const {
    double x = 0.0;
    double y = -(a_*x + c_)/b_;
    if (b_ == 0.0) {
        x = -c_/a_;
        y = 0.0;
    }
    return Point2(x, y);
}

Vector2 Line2::direction() const {
    return Vector2(b_, -a_);
}

Vector2 Line2::normal() const {
    return Vector2(a_, b_);
}

Line2 Line2::opposite() const {
    return Line2(-a_, -b_, -c_);
}

int Line2::side(const Point2& p) const {
    int result = 0;
    double x = p.getX();
    double y = p.getY();
    double dist = a_*x + b_*y + c_;
    if (dist > 0.0) {
        result = 1;
    } else if (dist < 0.0) {
        result = -1;
    }
    return result;
}

bool Line2::operator==(const Line2& l) const {
    bool result = false;
    if (a_ == l.a_ && b_ == l.b_ && c_ == l.c_) {
        result = true;
    } else if (normal().normalize() == l.normal().normalize() &&
            getC()/normal().length() == l.getC()/l.normal().length()) {
        result = true;
    }
    return result;
}

}
