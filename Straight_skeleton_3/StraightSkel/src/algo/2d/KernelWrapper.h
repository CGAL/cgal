/**
 * @file   algo/2d/KernelWrapper.h
 * @author Gernot Walzl
 * @date   2012-02-08
 */

#ifndef ALGO_2D_KERNELWRAPPER_H
#define ALGO_2D_KERNELWRAPPER_H

#include "config.h"
#ifdef USE_CGAL
    #include "cgal_kernel.h"
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_bisector.html
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_intersection.html
    #include <CGAL/intersections.h>
    // http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Kernel_23_ref/Function_squared_distance.html
    #include <CGAL/squared_distance_2.h>
#else
    #include "kernel/intersection.h"
    #include "kernel/distance.h"
    #include "kernel/bisector.h"
    #include "kernel/projection.h"
#endif

#include "debug.h"
#include "data/2d/ptrs.h"
#include "data/2d/KernelFactory.h"

namespace algo { namespace _2d {

using namespace data::_2d;

class KernelWrapper {
public:
    virtual ~KernelWrapper();

    static Point2SPtr intersection(Line2SPtr line1, Line2SPtr line2);
    static Line2SPtr bisector(Line2SPtr line1, Line2SPtr line2);
    static double distance(Point2SPtr p, Point2SPtr q);
    static double distance(Line2SPtr line, Point2SPtr point);
    static Line2SPtr opposite(Line2SPtr line);
    static Vector2SPtr perpendicular(Vector2SPtr vector);
    static Vector2SPtr normalize(Vector2SPtr vector);
    static Line2SPtr offsetLine(Line2SPtr line, double offset);
    static Point2SPtr offsetPoint(Point2SPtr point, Vector2SPtr dir, double offset);
    static double angle(Vector2SPtr v1, Vector2SPtr v2);
    static int side(Line2SPtr line, Point2SPtr point);

    static Point2SPtr projection(Line2SPtr line, Point2SPtr point);
    static int compatePoints(Vector2SPtr v_dir, Point2SPtr p_1, Point2SPtr p_2);
    static bool isInside(Point2SPtr p, Point2SPtr p_box_1, Point2SPtr p_box_2);

protected:
    KernelWrapper();
};

} }

#endif /* ALGO_2D_KERNELWRAPPER_H */

