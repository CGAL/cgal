/**
 * @file   kernel/projection.cpp
 * @author Gernot Walzl
 * @date   2013-01-29
 */

#include "kernel/projection.h"

#include <cmath>

namespace kernel {

Point2* projection(const Line2* line, const Point2* point) {
    Point2 p_line = line->point();
    Vector2 normal = line->normal().normalize();
    Point2 result = *point - (normal * ((*point - p_line) * normal));
    return new Point2(result);
}

Point3* projection(const Plane3* plane, const Point3* point) {
    Point3 p_plane = plane->point();
    Vector3 normal = plane->normal().normalize();
    Point3 result = *point - (normal * ((*point - p_plane) * normal));
    return new Point3(result);
}

Point3* projection(const Line3* line, const Point3* point) {
    Point3 p_line = line->point();
    Vector3 dir = line->direction().normalize();
    Point3 result = *point - ((*point - p_line) - (dir * ((*point - p_line) * dir)));
    return new Point3(result);
}

}
