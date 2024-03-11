/**
 * @file   algo/2d/KernelWrapper.cpp
 * @author Gernot Walzl
 * @date   2012-02-08
 */

#include "algo/2d/KernelWrapper.h"

#include <cmath>

namespace algo { namespace _2d {

KernelWrapper::KernelWrapper() {
    // intentionally does nothing
}

KernelWrapper::~KernelWrapper() {
    // intentionally does nothing
}

Point2SPtr KernelWrapper::intersection(Line2SPtr line1, Line2SPtr line2) {
    Point2SPtr result = Point2SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*line1, *line2);
    if (const CGAL::Point2 *ipoint = CGAL::object_cast<CGAL::Point2>(&obj)) {
        result = KernelFactory::createPoint2(*ipoint);
    }
#else
    result = Point2SPtr(kernel::intersection(&(*line1), &(*line2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Line2SPtr KernelWrapper::bisector(Line2SPtr line1, Line2SPtr line2) {
    Line2SPtr result = Line2SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createLine2(CGAL::bisector(*line1, *line2));
#else
    result = Line2SPtr(kernel::bisector(&(*line1), &(*line2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

double KernelWrapper::distance(Point2SPtr p, Point2SPtr q) {
    double result = 0.0;
#ifdef USE_CGAL
    result = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(*p, *q)));
#else
    result = kernel::distance(&(*p), &(*q));
#endif
    return result;
}

double KernelWrapper::distance(Line2SPtr line, Point2SPtr point) {
    double result = 0.0;
#ifdef USE_CGAL
    result = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(*line, *point)));
#else
    result = kernel::distance(&(*line), &(*point));
#endif
    return result;
}

Line2SPtr KernelWrapper::opposite(Line2SPtr line) {
    Line2SPtr result = KernelFactory::createLine2(line->opposite());
    DEBUG_SPTR(result);
    return result;
}

Vector2SPtr KernelWrapper::perpendicular(Vector2SPtr vector) {
    Vector2SPtr result;
    Vector2 v_perpend(-(*vector)[1], (*vector)[0]);
    result = KernelFactory::createVector2(v_perpend);
    return result;
}

Vector2SPtr KernelWrapper::normalize(Vector2SPtr vector) {
    Vector2SPtr result;
#ifdef USE_CGAL
    result = KernelFactory::createVector2(
            *vector / CGAL::sqrt(vector->squared_length()));
#else
    result = KernelFactory::createVector2(vector->normalize());
#endif
    return result;
}

Line2SPtr KernelWrapper::offsetLine(Line2SPtr line, double offset) {
    Line2SPtr result;
    Point2 p = line->point();
#ifdef USE_CGAL
    Vector2 v_dir = line->to_vector();
    Vector2 v_norm(-v_dir[1], v_dir[0]);
    Vector2 v_normal = v_norm / CGAL::sqrt(v_norm.squared_length());
#else
    Vector2 v_dir = line->direction();
    Vector2 v_normal = line->normal().normalize();
#endif
    Point2 p_trans = p + (v_normal * offset);
    Line2 line_trans(p_trans, v_dir);
    result = KernelFactory::createLine2(line_trans);
    DEBUG_SPTR(result);
    return result;
}

Point2SPtr KernelWrapper::offsetPoint(Point2SPtr point, Vector2SPtr dir, double offset) {
    Point2SPtr result;
#ifdef USE_CGAL
    Vector2 dir_normalized = *dir / CGAL::sqrt(dir->squared_length());
    Point2 p_moved = *point + (dir_normalized * offset);
#else
    Point2 p_moved = *point + (dir->normalize() * offset);
#endif
    result = KernelFactory::createPoint2(p_moved);
    DEBUG_SPTR(result);
    return result;
}

double KernelWrapper::angle(Vector2SPtr v1, Vector2SPtr v2) {
    double result = 0.0;
    double arg = 0.0;
#ifdef USE_CGAL
    arg = ((*v1)*(*v2)) / CGAL::sqrt(v1->squared_length() * v2->squared_length());
#else
    arg = ((*v1)*(*v2)) / sqrt(v1->squared_length() * v2->squared_length());
#endif
    // fixes issues with floating point precision
    if (arg <= -1.0) {
        result = M_PI;
    } else if (arg >= 1.0) {
        result = 0.0;
    } else {
        result = acos(arg);
    }
    return result;
}

int KernelWrapper::side(Line2SPtr line, Point2SPtr point) {
    int result = 0;
#ifdef USE_CGAL
    CGAL::Oriented_side side = line->oriented_side(*point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
#else
    result = line->side(*point);
#endif
    return result;
}


Point2SPtr KernelWrapper::projection(Line2SPtr line, Point2SPtr point) {
    Point2SPtr result = Point2SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPoint2(line->projection(*point));
#else
    result = Point2SPtr(kernel::projection(&(*line), &(*point)));
#endif
    DEBUG_SPTR(result);
    return result;
}

int KernelWrapper::compatePoints(Vector2SPtr v_dir, Point2SPtr p_1, Point2SPtr p_2) {
    int result = 0;
    double value = *v_dir * (*p_2 - *p_1);
    if (value > 0.0) {         // angle < M_PI/2.0
        result = -1;
    } else if (value < 0.0) {  // angle > M_PI/2.0
        result = 1;
    }
    return result;
}

bool KernelWrapper::isInside(Point2SPtr p, Point2SPtr p_box_1, Point2SPtr p_box_2) {
    bool result = true;
    for (unsigned int i = 0; i < 2; i++) {
        if ((*p_box_1)[i] < (*p_box_2)[i]) {
            if (!( (*p_box_1)[i] <= (*p)[i] && (*p)[i] <= (*p_box_2)[i] )) {
                result = false;
            }
        } else {
            if (!( (*p_box_2)[i] <= (*p)[i] && (*p)[i] <= (*p_box_1)[i] )) {
                result = false;
            }
        }
    }
    return result;
}

} }
