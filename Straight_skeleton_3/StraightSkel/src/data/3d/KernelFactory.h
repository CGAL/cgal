/**
 * @file   data/3d/KernelFactory.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_KERNELFACTORY_H
#define DATA_3D_KERNELFACTORY_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

class KernelFactory {
public:
    virtual ~KernelFactory();

    static Point3SPtr createPoint3(double x, double y, double z);
    static Point3SPtr createPoint3(const Point3& point);
    static Point3SPtr createPoint3(Vector3SPtr vector);

    /**
     * Returns the center of the sphere.
     */
    static Point3SPtr createPoint3(Sphere3SPtr sphere);

    static Segment3SPtr createSegment3(Point3SPtr src, Point3SPtr dst);
    static Segment3SPtr createSegment3(const Segment3& seg);

    static Line3SPtr createLine3(Point3SPtr p, Point3SPtr q);
    static Line3SPtr createLine3(const Line3& line);
    static Line3SPtr createLine3(Point3SPtr p, Vector3SPtr direction);

    static Plane3SPtr createPlane3(double a, double b, double c, double d);
    static Plane3SPtr createPlane3(Point3SPtr p, Point3SPtr q, Point3SPtr r);
    static Plane3SPtr createPlane3(Point3SPtr p, Vector3SPtr normal);
    static Plane3SPtr createPlane3(const Plane3& plane);

    static Sphere3SPtr createSphere3(Point3SPtr center, double radius);
    static Sphere3SPtr createSphere3(const Sphere3& sphere);

    static Vector3SPtr createVector3(double x, double y, double z);
    static Vector3SPtr createVector3(const Vector3& vector);
    static Vector3SPtr createVector3(Point3SPtr point);

    /**
     * Creates a vector with the same direction as the given line.
     */
    static Vector3SPtr createVector3(Line3SPtr line);

    /**
     * Returns the normal vector of the plane.
     */
    static Vector3SPtr createVector3(Plane3SPtr plane);

protected:
    KernelFactory();
};

} }

#endif /* DATA_3D_KERNELFACTORY_H */
