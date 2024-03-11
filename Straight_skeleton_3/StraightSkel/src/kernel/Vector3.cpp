/**
 * @file   kernel/Vector3.cpp
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#include "kernel/Vector3.h"

#include <cmath>
#include <stdexcept>

namespace kernel {

Vector3::Vector3() {
    for (int i = 0; i < 3; i++) {
        this->v_[i] = 0.0;
    }
}

Vector3::Vector3(double x, double y, double z) {
    this->v_[0] = x;
    this->v_[1] = y;
    this->v_[2] = z;
}

Vector3::Vector3(const Vector3& orig) {
    for (int i = 0; i < 3; i++) {
        this->v_[i] = orig.v_[i];
    }
}

Vector3::~Vector3() {
    // intentionally does nothing
}

double Vector3::operator[](unsigned int i) const {
    if (i >= 3) {
        throw std::out_of_range("Index out of bounds.");
    }
    return v_[i];
}

double Vector3::squared_length() const {
    double result = 0.0;
    for (int i = 0; i < 3; i++) {
        result += v_[i]*v_[i];
    }
    return result;
}

double Vector3::length() const {
    double result = sqrt(squared_length());
    return result;
}

Vector3 Vector3::normalize(void) const {
    double length = this->length();
    return Vector3(v_[0]/length, v_[1]/length, v_[2]/length);
}

double Vector3::angle(const Vector3& v) const {
    double result = acos( ((*this) * v) /
            sqrt(this->squared_length() * v.squared_length()) );
    return result;
}

Vector3 Vector3::operator+(const Vector3& v) const {
    return Vector3(v_[0]+v.v_[0], v_[1]+v.v_[1], v_[2]+v.v_[2]);
}

Vector3 Vector3::operator-(const Vector3& v) const {
    return Vector3(v_[0]-v.v_[0], v_[1]-v.v_[1], v_[2]-v.v_[2]);
}

Vector3 Vector3::operator*(double s) const {
    return Vector3(v_[0]*s, v_[1]*s, v_[2]*s);
}

Vector3 Vector3::operator/(double s) const {
    return Vector3(v_[0]/s, v_[1]/s, v_[2]/s);
}

double Vector3::operator*(const Vector3& v) const {
    double result = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        result += v_[i] * v.v_[i];
    }
    return result;
}

Vector3 Vector3::cross(const Vector3& v) const {
    double result_1 = v_[1]*v.v_[2] - v_[2]*v.v_[1];
    double result_2 = v_[2]*v.v_[0] - v_[0]*v.v_[2];
    double result_3 = v_[0]*v.v_[1] - v_[1]*v.v_[0];
    return Vector3(result_1, result_2, result_3);
}

bool Vector3::operator==(const Vector3& v) const {
    bool result = true;
    for (unsigned int i = 0; i < 3; i++) {
        if (v_[i] != v.v_[i]) {
            result = false;
            break;
        }
    }
    return result;
}

}
