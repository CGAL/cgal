/**
 * @file   kernel/bisector.cpp
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#include "kernel/bisector.h"

#include <cmath>

namespace kernel {

Line2* bisector(const Line2* line1, const Line2* line2) {
    Line2* result = 0;
    Point2* p = intersection(line1, line2);
    if (p) {
        Vector2 dir = line1->direction().normalize() +
                      line2->direction().normalize();
        result = new Line2(*p, dir);
        delete p;
    } else {
        Vector2 normal1 = line1->normal();
        Vector2 normal2 = line2->normal();
        if (normal1.angle(normal2) < M_PI/2.0) {
            unsigned int max_i = 0;
            double max_abs = 0.0;
            for (unsigned int i = 0; i < 2; i++) {
                double val_abs = normal2[i];
                if (val_abs < 0.0) {
                    val_abs *= -1.0;
                }
                if (val_abs >= max_abs) {
                    max_abs = val_abs;
                    max_i = i;
                }
            }
            double c = (line1->getC() +
                (line2->getC() * normal1[max_i]/normal2[max_i])) / 2.0;
            result = new Line2(normal1[0], normal1[1], c);
        }
    }
    return result;
}


Plane3* bisector(const Plane3* plane1, const Plane3* plane2) {
    Plane3* result = 0;
    Line3* l = intersection(plane1, plane2);
    if (l) {
        Vector3 normal = plane1->normal().normalize() +
                         plane2->normal().normalize();
        result = new Plane3(l->point(), normal);
//        Vector3 dir1 = (plane1->normal().cross(l->direction())).normalize();
//        Vector3 dir2 = (plane2->normal().cross(l->direction())).normalize();
//        Point3 p1 = l->point();
//        Point3 p2 = l->point() + l->direction();
//        Point3 p3 = l->point() + dir1 + dir2;
//        result = new Plane3(p1, p2, p3);
        delete l;
    } else {
        Vector3 normal1 = plane1->normal();
        Vector3 normal2 = plane2->normal();
        if (normal1.angle(normal2) < M_PI/2.0) {
            unsigned int max_i = 0;
            double max_abs = 0.0;
            for (unsigned int i = 0; i < 3; i++) {
                double val_abs = normal2[i];
                if (val_abs < 0.0) {
                    val_abs *= -1.0;
                }
                if (val_abs >= max_abs) {
                    max_abs = val_abs;
                    max_i = i;
                }
            }
            double d = (plane1->getD() +
                    (plane2->getD() * normal1[max_i]/normal2[max_i])) / 2.0;
            result = new Plane3(normal1[0], normal1[1], normal1[2], d);
        }
    }
    return result;
}

}
