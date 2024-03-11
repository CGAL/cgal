/**
 * @file   kernel/Line3.h
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#ifndef LINE3_H
#define LINE3_H

#include "kernel/Point3.h"
#include "kernel/Vector3.h"

namespace kernel {

class Line3 {
public:
    Line3();
    Line3(const Line3& orig);
    Line3(const Point3& p, const Point3& q);
    Line3(const Point3& p, const Vector3& dir);
    virtual ~Line3();
    Point3 point() const;
    Vector3 direction() const;
    Line3 opposite() const;
    bool hasOn(const Point3& point) const;
    bool operator==(const Line3& l) const;
protected:
    Point3* p_;
    Vector3* dir_;
};

}

#endif /* LINE3_H */

