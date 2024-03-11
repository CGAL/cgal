/**
 * @file   kernel/intersection.h
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Line3.h"
#include "kernel/Plane3.h"

namespace kernel {

Point2* intersection(const Line2* l1, const Line2* l2);

Point3* intersection(const Plane3* p1, const Plane3* p2, const Plane3* p3);
Line3* intersection(const Plane3* p1, const Plane3* p2);
Point3* intersection(const Plane3* plane, const Line3* line);

}

#endif /* INTERSECTION_H */

