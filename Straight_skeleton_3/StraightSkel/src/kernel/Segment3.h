/**
 * @file   kernel/Segment3.h
 * @author Gernot Walzl
 * @date   2011-11-11
 */

#ifndef SEGMENT3_H
#define SEGMENT3_H

#include "kernel/Point3.h"
#include "kernel/Line3.h"

namespace kernel {

class Segment3 {
public:
    Segment3(const Point3& p, const Point3& q);
    Segment3(const Segment3& orig);
    virtual ~Segment3();
    const Point3& getP() const;
    const Point3& getQ() const;
    void setP(const Point3& p);
    void setQ(const Point3& q);
    Line3 line() const;
protected:
    const Point3* p_;
    const Point3* q_;
};

}

#endif /* SEGMENT3_H */

