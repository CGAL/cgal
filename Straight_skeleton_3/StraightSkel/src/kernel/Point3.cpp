/**
 * @file   kernel/Point3.cpp
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#include "kernel/Point3.h"

#include <stdexcept>

namespace kernel {

Point3::Point3() {
    this->x_ = 0.0;
    this->y_ = 0.0;
    this->z_ = 0.0;
}

Point3::Point3(double x, double y, double z) {
    this->x_ = x;
    this->y_ = y;
    this->z_ = z;
}

Point3::Point3(const Point3& orig) {
    this->x_ = orig.x_;
    this->y_ = orig.y_;
    this->z_ = orig.z_;
}

Point3::~Point3() {
    // intentionally does nothing
}

double Point3::getX(void) const {
    return this->x_;
}

double Point3::getY(void) const {
    return this->y_;
}

double Point3::getZ(void) const {
    return this->z_;
}

double Point3::operator[](unsigned int i) const {
    if (i == 0) {
        return this->x_;
    } else if (i == 1) {
        return this->y_;
    } else if (i == 2) {
        return this->z_;
    } else {
        throw std::out_of_range("Index out of bounds.");
    }
}

Vector3 Point3::operator-(const Point3& q) const {
    double v_x = this->x_ - q.x_;
    double v_y = this->y_ - q.y_;
    double v_z = this->z_ - q.z_;
    return Vector3(v_x, v_y, v_z);
}

Point3 Point3::operator+(const Vector3& v) const {
    double x = this->x_ + v[0];
    double y = this->y_ + v[1];
    double z = this->z_ + v[2];
    return Point3(x, y, z);
}

Point3 Point3::operator-(const Vector3& v) const {
    double x = this->x_ - v[0];
    double y = this->y_ - v[1];
    double z = this->z_ - v[2];
    return Point3(x, y, z);
}

bool Point3::operator==(const Point3& q) const {
    return (this->x_ == q.getX() &&
            this->y_ == q.getY() &&
            this->z_ == q.getZ());
}

}
