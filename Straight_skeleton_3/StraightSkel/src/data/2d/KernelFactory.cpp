/**
 * @file   data/2d/KernelFactory.cpp
 * @author Gernot Walzl
 * @date   2011-11-25
 */

#include "data/2d/KernelFactory.h"

namespace data { namespace _2d {

KernelFactory::KernelFactory() {
    // intentionally does nothing
}

KernelFactory::~KernelFactory() {
    // intentionally does nothing
}

Point2SPtr KernelFactory::createPoint2(double x, double y) {
    return Point2SPtr(new Point2(x, y));
}

Point2SPtr KernelFactory::createPoint2(const Point2& point) {
    return Point2SPtr(new Point2(point));
}

Segment2SPtr KernelFactory::createSegment2(Point2SPtr src, Point2SPtr dst) {
    Segment2SPtr result = Segment2SPtr();
    if (src != dst) {
        result = Segment2SPtr(new Segment2(*src, *dst));
    }
    DEBUG_SPTR(result);
    return result;
}

Segment2SPtr KernelFactory::createSegment2(const Segment2& seg) {
    return Segment2SPtr(new Segment2(seg));
}

Line2SPtr KernelFactory::createLine2(Point2SPtr p, Point2SPtr q) {
    Line2SPtr result = Line2SPtr();
    if (p != q) {
        result = Line2SPtr(new Line2(*p, *q));
    }
    DEBUG_SPTR(result);
    return result;
}

Line2SPtr KernelFactory::createLine2(const Line2& line) {
    return Line2SPtr(new Line2(line));
}

Line2SPtr KernelFactory::createLine2(Segment2SPtr seg) {
    Line2SPtr result = Line2SPtr();
#ifdef USE_CGAL
    result = createLine2(seg->supporting_line());
#else
    result = createLine2(seg->line());
#endif
    return result;
}

Line2SPtr KernelFactory::createLine2(Point2SPtr p, Vector2SPtr direction) {
    return Line2SPtr(new Line2(*p, *direction));
}

Vector2SPtr KernelFactory::createVector2(double x, double y) {
    return Vector2SPtr(new Vector2(x, y));
}

Vector2SPtr KernelFactory::createVector2(const Vector2& vector) {
    return Vector2SPtr(new Vector2(vector));
}

Vector2SPtr KernelFactory::createVector2(Point2SPtr point) {
    Vector2SPtr result;
#ifdef USE_CGAL
    result = Vector2SPtr(new Vector2(point->x(), point->y()));
#else
    result = Vector2SPtr(new Vector2(point->getX(), point->getY()));
#endif
    return result;
}

Vector2SPtr KernelFactory::createVector2(Line2SPtr line) {
    Vector2SPtr result;
#ifdef USE_CGAL
    result = createVector2(line->to_vector());
#else
    result = createVector2(line->direction());
#endif
    return result;
}

} }
