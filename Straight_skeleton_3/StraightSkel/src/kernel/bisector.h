/**
 * @file   kernel/bisector.h
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#ifndef BISECTOR_H
#define BISECTOR_H

#include "kernel/Point2.h"
#include "kernel/Line2.h"
#include "kernel/Point3.h"
#include "kernel/Plane3.h"
#include "kernel/intersection.h"

namespace kernel {

Line2* bisector(const Line2* line1, const Line2* line2);

Plane3* bisector(const Plane3* plane1, const Plane3* plane2);

}

#endif /* BISECTOR_H */

