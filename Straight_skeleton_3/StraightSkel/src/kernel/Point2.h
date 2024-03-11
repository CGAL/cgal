/**
 * @file   kernel/Point2.h
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#ifndef POINT2_H
#define POINT2_H

#include "kernel/Vector2.h"

namespace kernel {

class Point2 {
public:
    Point2();
    Point2(double x, double y);
    Point2(const Point2& orig);
    virtual ~Point2();
    double getX(void) const;
    double getY(void) const;
    double operator[](unsigned int i) const;
    Vector2 operator-(const Point2& q) const;
    Point2 operator+(const Vector2& v) const;
    Point2 operator-(const Vector2& v) const;
    bool operator==(const Point2& q) const;
protected:
    double x_;
    double y_;
};

}

#endif /* POINT2_H */
