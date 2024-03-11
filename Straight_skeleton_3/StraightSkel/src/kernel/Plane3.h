/**
 * @file   kernel/Plane3.h
 * @author Gernot Walzl
 * @date   2011-11-14
 */

#ifndef PLANE3_H
#define PLANE3_H

#include "kernel/Point3.h"
#include "kernel/Vector3.h"

namespace kernel {

/*!
 * a*x + b*y + c*z + d = 0
 */
class Plane3 {
public:
    Plane3();
    Plane3(const Plane3& orig);
    Plane3(const Point3& p, const Point3& q, const Point3& r);
    Plane3(const Point3& p, const Vector3& normal);
    Plane3(double a, double b, double c, double d);
    virtual ~Plane3();
    double getA() const;
    double getB() const;
    double getC() const;
    double getD() const;
    Vector3 normal() const;
    Point3 point() const;
    Plane3 opposite() const;
    int side(const Point3& p) const;
    bool operator==(const Plane3& plane) const;
protected:
    double a_;
    double b_;
    double c_;
    double d_;
};

}

#endif /* PLANE3_H */

