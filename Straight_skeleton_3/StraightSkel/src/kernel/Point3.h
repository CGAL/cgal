/**
 * @file   kernel/Point3.h
 * @author Gernot Walzl
 * @date   2011-11-10
 */

#ifndef POINT3_H
#define POINT3_H

#include "kernel/Vector3.h"

namespace kernel {

class Point3 {
public:
    Point3();
    Point3(double x, double y, double z);
    Point3(const Point3& orig);
    virtual ~Point3();
    double getX(void) const;
    double getY(void) const;
    double getZ(void) const;
    double operator[](unsigned int i) const;
    Vector3 operator-(const Point3& q) const;
    Point3 operator+(const Vector3& v) const;
    Point3 operator-(const Vector3& v) const;
    bool operator==(const Point3& q) const;
protected:
    double x_;
    double y_;
    double z_;
};

}

#endif /* POINT3_H */
