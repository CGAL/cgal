/**
 * @file   kernel/Line3.cpp
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#include "kernel/Line3.h"

namespace kernel {

Line3::Line3() {
    this->p_ = new Point3();
    this->dir_ = new Vector3();
}

Line3::Line3(const Line3& orig) {
    this->p_ = new Point3(orig.point());
    this->dir_ = new Vector3(orig.direction());
}

Line3::Line3(const Point3& p, const Point3& q) {
    this->p_ = new Point3(p);
    this->dir_ = new Vector3(q-p);
}

Line3::Line3(const Point3& p, const Vector3& dir) {
    this->p_ = new Point3(p);
    this->dir_ = new Vector3(dir);
}

Line3::~Line3() {
    delete this->p_;
    delete this->dir_;
}

Point3 Line3::point() const {
    return Point3(*p_);
}

Vector3 Line3::direction() const {
    return Vector3(*dir_);
}

Line3 Line3::opposite() const {
    return Line3(*p_, (*dir_)*(-1.0));
}

bool Line3::hasOn(const Point3& point) const {
    Point3 p1 = *p_;
    Point3 p2 = *p_ + *dir_;
    double distance = ((point-p1).cross(point-p2)).length() / (p2-p1).length();
    return (distance == 0.0);
}

bool Line3::operator==(const Line3& l) const {
    bool result = (direction().normalize() == l.direction().normalize() &&
            hasOn(l.point()));
    return result;
}

}
