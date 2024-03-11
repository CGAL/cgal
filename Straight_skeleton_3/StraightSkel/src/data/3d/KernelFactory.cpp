/**
 * @file   data/3d/KernelFactory.cpp
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#include "data/3d/KernelFactory.h"

#include "debug.h"

namespace data { namespace _3d {

KernelFactory::KernelFactory() {
    // intentionally does nothing
}

KernelFactory::~KernelFactory() {
    // intentionally does nothing
}

Point3SPtr KernelFactory::createPoint3(double x, double y, double z) {
    return Point3SPtr(new Point3(x, y, z));
}

Point3SPtr KernelFactory::createPoint3(const Point3& point) {
    return Point3SPtr(new Point3(point));
}

Point3SPtr KernelFactory::createPoint3(Vector3SPtr vector) {
    Point3SPtr result;
    result = Point3SPtr(new Point3((*vector)[0], (*vector)[1], (*vector)[2]));
    return result;
}

Point3SPtr KernelFactory::createPoint3(Sphere3SPtr sphere) {
    Point3SPtr result;
#ifdef USE_CGAL
    result = createPoint3(sphere->center());
#else
    result = createPoint3(sphere->getCenter());
#endif
    return result;
}

Segment3SPtr KernelFactory::createSegment3(Point3SPtr src, Point3SPtr dst) {
    Segment3SPtr result = Segment3SPtr();
    if (src != dst) {
        result = Segment3SPtr(new Segment3(*src, *dst));
    }
    DEBUG_SPTR(result)
    return result;
}

Segment3SPtr KernelFactory::createSegment3(const Segment3& seg) {
    return Segment3SPtr(new Segment3(seg));
}

Line3SPtr KernelFactory::createLine3(Point3SPtr p, Point3SPtr q) {
    Line3SPtr result = Line3SPtr();
    if (p != q) {
        result = Line3SPtr(new Line3(*p, *q));
    }
    DEBUG_SPTR(result);
    return result;
}

Line3SPtr KernelFactory::createLine3(const Line3& line) {
    return Line3SPtr(new Line3(line));
}

Line3SPtr KernelFactory::createLine3(Point3SPtr p, Vector3SPtr direction) {
    return Line3SPtr(new Line3(*p, *direction));
}

Plane3SPtr KernelFactory::createPlane3(double a, double b, double c, double d) {
    return Plane3SPtr(new Plane3(a, b, c, d));
}

Plane3SPtr KernelFactory::createPlane3(Point3SPtr p, Point3SPtr q, Point3SPtr r) {
    Plane3SPtr result = Plane3SPtr();
    if (p != q && q != r && r != p) {
        result = Plane3SPtr(new Plane3(*p, *q, *r));
    }
    DEBUG_SPTR(result);
    return result;
}

Plane3SPtr KernelFactory::createPlane3(Point3SPtr p, Vector3SPtr normal) {
    return Plane3SPtr(new Plane3(*p, *normal));
}

Plane3SPtr KernelFactory::createPlane3(const Plane3& plane) {
    return Plane3SPtr(new Plane3(plane));
}

Sphere3SPtr KernelFactory::createSphere3(Point3SPtr center, double radius) {
    return Sphere3SPtr(new Sphere3(*center, radius));
}

Sphere3SPtr KernelFactory::createSphere3(const Sphere3& sphere) {
    return Sphere3SPtr(new Sphere3(sphere));
}

Vector3SPtr KernelFactory::createVector3(double x, double y, double z) {
    return Vector3SPtr(new Vector3(x, y, z));
}

Vector3SPtr KernelFactory::createVector3(const Vector3& vector) {
    return Vector3SPtr(new Vector3(vector));
}

Vector3SPtr KernelFactory::createVector3(Point3SPtr point) {
    Vector3SPtr result;
#ifdef USE_CGAL
    result = Vector3SPtr(new Vector3(
            point->x(), point->y(), point->z()));
#else
    result = Vector3SPtr(new Vector3(
            point->getX(), point->getY(), point->getZ()));
#endif
    return result;
}

Vector3SPtr KernelFactory::createVector3(Line3SPtr line) {
    Vector3SPtr result;
#ifdef USE_CGAL
    result = createVector3(line->to_vector());
#else
    result = createVector3(line->direction());
#endif
    return result;
}

Vector3SPtr KernelFactory::createVector3(Plane3SPtr plane) {
    Vector3SPtr result;
#ifdef USE_CGAL
    result = createVector3(plane->orthogonal_vector());
#else
    result = createVector3(plane->normal());
#endif
    return result;
}

} }
