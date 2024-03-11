/**
 * @file   kernel/intersection.cpp
 * @author Gernot Walzl
 * @date   2012-02-25
 */

#include "kernel/intersection.h"

#include <cmath>

namespace kernel {

Point2* intersection(const Line2* l1, const Line2* l2) {
    if (l1 == l2) {
        return 0;
    }
    if (*l1 == *l2) {
        return 0;
    }
    double det = l1->getA() * l2->getB() - l2->getA() * l1->getB();
    if (det == 0.0) {
        return 0;
    }
    double p_x = (l1->getB() * l2->getC() - l1->getC() * l2->getB()) / det;
    double p_y = (l1->getC() * l2->getA() - l2->getC() * l1->getA()) / det;
    Point2* result = new Point2(p_x, p_y);
    return result;
}


Point3* intersection(const Plane3* p1, const Plane3* p2, const Plane3* p3) {
    if (p1 == p2 || p2 == p3 || p3 == p1) {
        return 0;
    }
    if (*p1 == *p2 || *p2 == *p3 || *p3 == *p1) {
        return 0;
    }
    // Sarrus
    double det = p1->getA() * p2->getB() * p3->getC()
            + p1->getB() * p2->getC() * p3->getA()
            + p1->getC() * p2->getA() * p3->getB()
            - p3->getA() * p2->getB() * p1->getC()
            - p3->getB() * p2->getC() * p1->getA()
            - p3->getC() * p2->getA() * p1->getB();
    if (det == 0.0) {
        return 0;
    }
    double p_x = (-(p2->getB()*p3->getC() - p2->getC()*p3->getB()) * p1->getD()
            - (p1->getC()*p3->getB() - p1->getB()*p3->getC()) * p2->getD()
            - (p1->getB()*p2->getC() - p1->getC()*p2->getB()) * p3->getD())
                / det;
    double p_y = (-(p2->getC()*p3->getA() - p2->getA()*p3->getC()) * p1->getD()
            - (p1->getA()*p3->getC() - p1->getC()*p3->getA()) * p2->getD()
            - (p1->getC()*p2->getA() - p1->getA()*p2->getC()) * p3->getD())
                / det;
    double p_z = (-(p2->getA()*p3->getB() - p2->getB()*p3->getA()) * p1->getD()
            - (p1->getB()*p3->getA() - p1->getA()*p3->getB()) * p2->getD()
            - (p1->getA()*p2->getB() - p1->getB()*p2->getA()) * p3->getD())
                / det;
    Point3* result = new Point3(p_x, p_y, p_z);
    return result;
}

Line3* intersection(const Plane3* p1, const Plane3* p2) {
    if (p1 == p2) {
        return 0;
    }
    if (*p1 == *p2) {
        return 0;
    }
    Vector3 v_dir = p1->normal().cross(p2->normal());
    if (v_dir.squared_length() == 0.0) {  // parallel
        return 0;
    }
    unsigned int max_i = 0;
    double max_abs = 0.0;
    for (unsigned int i = 0; i < 3; i++) {
        double val_abs = v_dir[i];
        if (val_abs < 0.0) {
            val_abs *= -1.0;
        }
        if (val_abs >= max_abs) {
            max_abs = val_abs;
            max_i = i;
        }
    }
    double det = 0.0;
    double p_x = 0.0;
    double p_y = 0.0;
    double p_z = 0.0;
    if (max_i == 2) {
        p_z = 0.0;
        det = p1->getA()*p2->getB() - p2->getA()*p1->getB();
        p_x = (p1->getB()*p2->getD() - p2->getB()*p1->getD()) / det;
        p_y = (p2->getA()*p1->getD() - p1->getA()*p2->getD()) / det;
    } else if (max_i == 1) {
        p_y = 0.0;
        det = p1->getA()*p2->getC() - p2->getA()*p1->getC();
        p_x = (p1->getC()*p2->getD() - p2->getC()*p1->getD()) / det;
        p_z = (p2->getA()*p1->getD() - p1->getA()*p2->getD()) / det;
    } else if (max_i == 0) {
        p_x = 0.0;
        det = p1->getB()*p2->getC() - p2->getB()*p1->getC();
        p_y = (p1->getC()*p2->getD() - p2->getC()*p1->getD()) / det;
        p_z = (p2->getB()*p1->getD() - p1->getB()*p2->getD()) / det;
    }
    Point3 p(p_x, p_y, p_z);
    Line3* result = new Line3(p, v_dir);
    return result;
}

Point3* intersection(const Plane3* plane, const Line3* line) {
    Point3 p0 = line->point();
    double skalprod = plane->normal() * line->direction();
    if (skalprod == 0.0) {
        return 0;
    }
    double lambda = -(plane->getA()*p0.getX() + plane->getB()*p0.getY() +
                      plane->getC()*p0.getZ() + plane->getD()) / skalprod;
    Point3* result = new Point3(p0 + line->direction()*lambda);
    return result;
}

}
