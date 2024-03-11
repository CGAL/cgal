/**
 * @file   kernel/Plane3.cpp
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#include "kernel/Plane3.h"

namespace kernel {

Plane3::Plane3() {
    this->a_ = 0.0;
    this->b_ = 0.0;
    this->c_ = 0.0;
    this->d_ = 0.0;
}

Plane3::Plane3(const Plane3& orig) {
    this->a_ = orig.a_;
    this->b_ = orig.b_;
    this->c_ = orig.c_;
    this->d_ = orig.d_;
}

/*!
 *                        ^(-1)
 *     [a]   [x1  y1  z1]       [-1]
 * 1/d [b] = [x2  y2  z2]       [-1]
 *     [c]   [x3  y3  z3]       [-1]
 */
Plane3::Plane3(const Point3& p1, const Point3& p2, const Point3& p3) {
    a_ = p2.getY()*p3.getZ() + p1.getZ()*p3.getY() + p1.getY()*p2.getZ()
            - p2.getZ()*p3.getY() - p1.getY()*p3.getZ() - p1.getZ()*p2.getY();
    b_ = p2.getZ()*p3.getX() + p1.getX()*p3.getZ() + p1.getZ()*p2.getX()
            - p2.getX()*p3.getZ() - p1.getZ()*p3.getX() - p1.getX()*p2.getZ();
    c_ = p2.getX()*p3.getY() + p1.getY()*p3.getX() + p1.getX()*p2.getY()
            - p2.getY()*p3.getX() - p1.getX()*p3.getY() - p1.getY()*p2.getX();
    d_ = - a_*p1.getX() - b_*p1.getY() - c_*p1.getZ();
}

Plane3::Plane3(const Point3& p, const Vector3& normal) {
    a_ = normal[0];
    b_ = normal[1];
    c_ = normal[2];
    d_ = - a_*p.getX() - b_*p.getY() - c_*p.getZ();
}

Plane3::Plane3(double a, double b, double c, double d) {
    this->a_ = a;
    this->b_ = b;
    this->c_ = c;
    this->d_ = d;
}

Plane3::~Plane3() {
    // intentionally does nothing
}

double Plane3::getA() const {
    return this->a_;
}

double Plane3::getB() const {
    return this->b_;
}

double Plane3::getC() const {
    return this->c_;
}

double Plane3::getD() const {
    return this->d_;
}

Vector3 Plane3::normal() const {
    return Vector3(a_, b_, c_);
}

Point3 Plane3::point() const {
    double x = 0.0;
    double y = 0.0;
    double z = -(a_*x + b_*y + d_) / c_;
    if (c_ == 0.0) {
        z = 0.0;
        y = -(a_*x + d_) / b_;
        if (b_ == 0.0) {
            x = -d_/a_;
            y = 0.0;
        }
    }
    return Point3(x, y, z);
}

Plane3 Plane3::opposite() const {
    return Plane3(-a_, -b_, -c_, -d_);
}

int Plane3::side(const Point3& p) const {
    int result = 0;
    double x = p.getX();
    double y = p.getY();
    double z = p.getZ();
    double dist = a_*x + b_*y + c_*z + d_;
    if (dist > 0.0) {
        result = 1;
    } else if (dist < 0.0) {
        result = -1;
    }
    return result;
}

bool Plane3::operator==(const Plane3& plane) const {
    bool result = false;
    if (a_ == plane.a_ && b_ == plane.b_ && c_ == plane.c_ && d_ == plane.d_) {
        result = true;
    } else if (normal().normalize() == plane.normal().normalize() &&
            getD()/normal().length() == plane.getD()/plane.normal().length()) {
        result = true;
    }
    return result;
}

}
