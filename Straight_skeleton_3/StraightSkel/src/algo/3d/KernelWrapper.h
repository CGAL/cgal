/**
 * @file   algo/3d/KernelWrapper.h
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#ifndef ALGO_3D_KERNELWRAPPER_H
#define ALGO_3D_KERNELWRAPPER_H

#include <cmath>

#include "config.h"
#ifdef USE_CGAL
    #include "cgal_kernel.h"
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_bisector.html
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_intersection.html
    #include <CGAL/intersections.h>
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_squared_distance.html
    #include <CGAL/squared_distance_3.h>
#else
    #include "kernel/intersection.h"
    #include "kernel/distance.h"
    #include "kernel/bisector.h"
    #include "kernel/projection.h"
#endif

#include "debug.h"
#include "data/3d/ptrs.h"
#include "data/3d/KernelFactory.h"

namespace algo { namespace _3d {

using namespace data::_3d;

class KernelWrapper {
public:
    virtual ~KernelWrapper();

    static Point3SPtr intersection(Plane3SPtr plane1, Plane3SPtr plane2, Plane3SPtr plane3);
    static Line3SPtr intersection(Plane3SPtr plane1, Plane3SPtr plane2);
    static Point3SPtr intersection(Plane3SPtr plane, Line3SPtr line);

    /**
     * If a line intersects a sphere, there are 2 intersection points.
     * The first one is returned here.
     */
    static Point3SPtr intersection(Sphere3SPtr sphere, Line3SPtr line);

    static Plane3SPtr bisector(Plane3SPtr plane1, Plane3SPtr plane2);

    static double distance(Point3SPtr p1, Point3SPtr p2);
    static double distance(Plane3SPtr plane, Point3SPtr point);
    static double distance(Line3SPtr line, Point3SPtr point);

    static Plane3SPtr opposite(Plane3SPtr plane);
    static Line3SPtr opposite(Line3SPtr line);

    static Vector3SPtr normalize(Vector3SPtr v);

    static Plane3SPtr offsetPlane(Plane3SPtr plane, double offset);
    static Point3SPtr offsetPoint(Point3SPtr point, Vector3SPtr dir, double offset);

    /**
     * http://de.wikipedia.org/wiki/Drehmatrix
     *
     *              [ n_x^2 (1 - cos(alpha)) + cos(alpha)         n_x n_y (1 - cos(alpha)) - n_z sin(alpha)   n_x n_z (1 - cos(alpha)) + n_y sin(alpha) ]
     * R_n(alpha) = [ n_y n_x (1 - cos(alpha)) + n_z sin(alpha)   n_y^2 (1 - cos(alpha)) + cos(alpha)         n_y n_z (1 - cos(alpha)) - n_x sin(alpha) ]
     *              [ n_z n_x (1 - cos(alpha)) - n_y sin(alpha)   n_z n_y (1 - cos(alpha)) + n_x sin(alpha)   n_z^2 (1 - cos(alpha)) + cos(alpha)       ]
     */
    static Vector3SPtr rotateVector(Vector3SPtr vector, Vector3SPtr axis, double angle);
    static Plane3SPtr rotatePlane(Plane3SPtr plane, Line3SPtr line, double angle);

    static int side(Plane3SPtr plane, Point3SPtr point);
    static int orientation(Line3SPtr line1, Line3SPtr line2);

    static double angle(Vector3SPtr v1, Vector3SPtr v2);
    static double angle(Line3SPtr line1, Line3SPtr line2);
    static double angle(Plane3SPtr plane, Line3SPtr line);

    /**
     * Computes the angle between the normal vectors of given planes.
     */
    static double angle(Plane3SPtr plane1, Plane3SPtr plane2);

    static bool isInside(Point3SPtr p, Point3SPtr p_box_1, Point3SPtr p_box_2);

    static Vector3SPtr cross(Vector3SPtr v1, Vector3SPtr v2);

    static Point3SPtr projection(Line3SPtr line, Point3SPtr point);
    static Point3SPtr projection(Plane3SPtr plane, Point3SPtr point);

    static int comparePoints(Vector3SPtr v_dir, Point3SPtr p_1, Point3SPtr p_2);

    static Point3SPtr replaceCoord(Point3SPtr point, Point3SPtr replacement,
            unsigned int coord);

protected:
    KernelWrapper();
};

} }

#endif /* ALGO_3D_KERNELWRAPPER_H */
