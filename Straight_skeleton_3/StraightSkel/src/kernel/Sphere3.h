/**
 * @file   kernel/Sphere3.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef SPHERE3_H
#define SPHERE3_H

#include "kernel/Point3.h"

namespace kernel {

class Sphere3 {
public:
    Sphere3(const Sphere3& orig);
    Sphere3(const Point3& center, double radius);
    virtual ~Sphere3();
    const Point3& getCenter() const;
    double getRadius() const;
    void setCenter(const Point3& center);
    void setRadius(double radius);
protected:
    const Point3* center_;
    double radius_;
};

}

#endif /* SPHERE3_H */

