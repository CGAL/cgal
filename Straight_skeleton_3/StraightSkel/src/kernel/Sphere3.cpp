/**
 * @file   kernel/Sphere3.cpp
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#include "kernel/Sphere3.h"

namespace kernel {

Sphere3::Sphere3(const Point3& center, double radius) {
    this->center_ = &center;
    this->radius_ = radius;
}

Sphere3::Sphere3(const Sphere3& orig) {
    this->center_ = orig.center_;
    this->radius_ = orig.radius_;
}

Sphere3::~Sphere3() {
    // intentionally does nothing
}

const Point3& Sphere3::getCenter() const {
    return *(this->center_);
}

double Sphere3::getRadius() const {
    return this->radius_;
}

void Sphere3::setCenter(const Point3& center) {
    this->center_ = &center;
}

void Sphere3::setRadius(double radius) {
    this->radius_ = radius;
}

}
