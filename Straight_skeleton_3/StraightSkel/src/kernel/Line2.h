/**
 * @file   kernel/Line2.h
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#ifndef LINE2_H
#define LINE2_H

#include "kernel/Point2.h"
#include "kernel/Vector2.h"

namespace kernel {

/*!
 * a*x + b*y + c = 0
 */
class Line2 {
public:
    Line2();
    Line2(const Line2& orig);
    Line2(const Point2& p, const Point2& q);
    Line2(const Point2& p, const Vector2& direction);
    Line2(double a, double b, double c);
    virtual ~Line2();
    double getA() const;
    double getB() const;
    double getC() const;
    Point2 point() const;
    Vector2 direction() const;
    Vector2 normal() const;
    Line2 opposite() const;
    int side(const Point2& p) const;
    bool operator==(const Line2& l) const;
protected:
    double a_;
    double b_;
    double c_;
};

}

#endif /* LINE2_H */

