/**
 * @file   kernel/Vector2.cpp
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#include "kernel/Vector2.h"

#include <cmath>
#include <stdexcept>

namespace kernel {

Vector2::Vector2() {
    for (int i = 0; i < 2; i++) {
        this->v_[i] = 0.0;
    }
}

Vector2::Vector2(double x, double y) {
    this->v_[0] = x;
    this->v_[1] = y;
}

Vector2::Vector2(const Vector2& orig) {
    for (int i = 0; i < 2; i++) {
        this->v_[i] = orig.v_[i];
    }
}

Vector2::~Vector2() {
    // intentionally does nothing
}

double Vector2::operator[](unsigned int i) const {
    if (i >= 2) {
        throw std::out_of_range("Index out of bounds.");
    }
    return v_[i];
}

double Vector2::squared_length(void) const {
    double result = 0.0;
    for (int i = 0; i < 2; i++) {
        result += v_[i]*v_[i];
    }
    return result;
}

double Vector2::length(void) const {
    double result = sqrt(squared_length());
    return result;
}

Vector2 Vector2::normalize(void) const {
    double length = this->length();
    return Vector2(v_[0]/length, v_[1]/length);
}

double Vector2::angle(const Vector2& v) const {
    double result = acos( ((*this) * v) /
            sqrt(this->squared_length() * v.squared_length()) );
    return result;
}

Vector2 Vector2::operator+(const Vector2& v) const {
    return Vector2(v_[0]+v.v_[0], v_[1]+v.v_[1]);
}

Vector2 Vector2::operator-(const Vector2& v) const {
    return Vector2(v_[0]-v.v_[0], v_[1]-v.v_[1]);
}

Vector2 Vector2::operator*(double s) const {
    return Vector2(v_[0]*s, v_[1]*s);
}

Vector2 Vector2::operator/(double s) const {
    return Vector2(v_[0]/s, v_[1]/s);
}

double Vector2::operator*(const Vector2& v) const {
    double result = 0.0;
    for (unsigned int i = 0; i < 2; i++) {
        result += v_[i] * v.v_[i];
    }
    return result;
}

bool Vector2::operator==(const Vector2& v) const {
    bool result = true;
    for (unsigned int i = 0; i < 2; i++) {
        if (v_[i] != v.v_[i]) {
            result = false;
            break;
        }
    }
    return result;
}

}
