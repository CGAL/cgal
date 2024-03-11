/**
 * @file   kernel/distance.cpp
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#include "kernel/distance.h"

#include <cmath>

namespace kernel {

double distance(const Point2* p, const Point2* q) {
    double dx = q->getX() - p->getX();
    double dy = q->getY() - p->getY();
    double result = sqrt(dx*dx + dy*dy);
    return result;
}

double distance(const Line2* line, const Point2* point) {
    double a = line->getA();
    double b = line->getB();
    double c = line->getC();
    double x = point->getX();
    double y = point->getY();
    double result = (a*x + b*y + c)/sqrt(a*a + b*b);
    if (result < 0.0) {
        result *= -1.0;
    }
    return result;
}


double distance(const Point3* p, const Point3* q) {
    double dx = q->getX() - p->getX();
    double dy = q->getY() - p->getY();
    double dz = q->getZ() - p->getZ();
    double result = sqrt(dx*dx + dy*dy + dz*dz);
    return result;
}

double distance(const Plane3* plane, const Point3* point) {
    double a = plane->getA();
    double b = plane->getB();
    double c = plane->getC();
    double d = plane->getD();
    double x = point->getX();
    double y = point->getY();
    double z = point->getZ();
    double result = (a*x + b*y + c*z + d)/sqrt(a*a + b*b + c*c);
    if (result < 0.0) {
        result *= -1.0;
    }
    return result;
}

double distance(const Line3* line, const Point3* point) {
    Point3 p1 = line->point();
    Point3 p2 = p1 + line->direction();
    double result = (((*point)-p1).cross((*point)-p2)).length() / (p2-p1).length();
    return result;
}

}
