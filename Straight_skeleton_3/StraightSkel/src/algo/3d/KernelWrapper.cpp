/**
 * @file   algo/3d/KernelWrapper.cpp
 * @author Gernot Walzl
 * @date   2012-03-08
 */

#include "algo/3d/KernelWrapper.h"

namespace algo { namespace _3d {

KernelWrapper::KernelWrapper() {
    // intentionally does nothing
}

KernelWrapper::~KernelWrapper() {
    // intentionally does nothing
}

Point3SPtr KernelWrapper::intersection(Plane3SPtr plane1, Plane3SPtr plane2, Plane3SPtr plane3) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane1, *plane2);
    if (const CGAL::Line3 *iline = CGAL::object_cast<CGAL::Line3>(&obj)) {
        CGAL::Object obj = CGAL::intersection(*iline, *plane3);
        if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
            result = KernelFactory::createPoint3(*ipoint);
        }
    }
#else
    result = Point3SPtr(kernel::intersection(&(*plane1), &(*plane2), &(*plane3)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Line3SPtr KernelWrapper::intersection(Plane3SPtr plane1, Plane3SPtr plane2) {
    Line3SPtr result = Line3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane1, *plane2);
    if (const CGAL::Line3 *iline = CGAL::object_cast<CGAL::Line3>(&obj)) {
        result = KernelFactory::createLine3(*iline);
    }
#else
    result = Line3SPtr(kernel::intersection(&(*plane1), &(*plane2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::intersection(Plane3SPtr plane, Line3SPtr line) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    CGAL::Object obj = CGAL::intersection(*plane, *line);
    if (const CGAL::Point3 *ipoint = CGAL::object_cast<CGAL::Point3>(&obj)) {
        result = KernelFactory::createPoint3(*ipoint);
    }
#else
    result = Point3SPtr(kernel::intersection(&(*plane), &(*line)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::intersection(Sphere3SPtr sphere, Line3SPtr line) {
    Point3SPtr result = Point3SPtr();
    Point3SPtr p_center = KernelFactory::createPoint3(sphere);
    double radius;
#ifdef USE_CGAL
    radius = CGAL::sqrt(sphere->squared_radius());
#else
    radius = sphere->getRadius();
#endif
    Vector3SPtr dir = normalize(KernelFactory::createVector3(line));
    Plane3SPtr plane = KernelFactory::createPlane3(p_center, dir);
    Point3SPtr p_intersect = intersection(plane, line);
    double dist = distance(p_center, p_intersect);
    if (dist == radius) {
        result = p_intersect;
    } else if (dist < radius) {
        double amount = -sqrt(radius*radius - dist*dist);
        result = KernelFactory::createPoint3((*p_intersect) + ((*dir)*amount));
    }
    //DEBUG_SPTR(result);
    return result;
}

Plane3SPtr KernelWrapper::bisector(Plane3SPtr plane1, Plane3SPtr plane2) {
    Plane3SPtr result = Plane3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPlane3(CGAL::bisector(*plane1, *plane2));
#else
    result = Plane3SPtr(kernel::bisector(&(*plane1), &(*plane2)));
#endif
    DEBUG_SPTR(result);
    return result;
}

double KernelWrapper::distance(Point3SPtr p1, Point3SPtr p2) {
    double result = 0.0;
#ifdef USE_CGAL
    result = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(*p1, *p2)));
#else
    result = kernel::distance(&(*p1), &(*p2));
#endif
    return result;
}

double KernelWrapper::distance(Plane3SPtr plane, Point3SPtr point) {
    double result = 0.0;
#ifdef USE_CGAL
    result = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(*plane, *point)));
#else
    result = kernel::distance(&(*plane), &(*point));
#endif
    return result;
}

double KernelWrapper::distance(Line3SPtr line, Point3SPtr point) {
    double result = 0.0;
#ifdef USE_CGAL
    result = CGAL::to_double(CGAL::sqrt(CGAL::squared_distance(*line, *point)));
#else
    result = kernel::distance(&(*line), &(*point));
#endif
    return result;
}

Plane3SPtr KernelWrapper::opposite(Plane3SPtr plane) {
    Plane3SPtr result = KernelFactory::createPlane3(plane->opposite());
    DEBUG_SPTR(result);
    return result;
}

Line3SPtr KernelWrapper::opposite(Line3SPtr line) {
    Line3SPtr result = KernelFactory::createLine3(line->opposite());
    DEBUG_SPTR(result);
    return result;
}

Vector3SPtr KernelWrapper::normalize(Vector3SPtr v) {
    Vector3SPtr result;
#ifdef USE_CGAL
    result = KernelFactory::createVector3(*v / CGAL::sqrt(v->squared_length()));
#else
    result = KernelFactory::createVector3(v->normalize());
#endif
    return result;
}

Plane3SPtr KernelWrapper::offsetPlane(Plane3SPtr plane, double offset) {
    Plane3SPtr result = Plane3SPtr();
    Point3 p = plane->point();
#ifdef USE_CGAL
    Vector3 v_norm = plane->orthogonal_vector();
    Vector3 v_normal = v_norm / CGAL::sqrt(v_norm.squared_length());
#else
    Vector3 v_normal = plane->normal().normalize();
#endif
    Point3 p_trans = p + (v_normal * offset);
    Plane3 plane_trans(p_trans, v_normal);
    result = KernelFactory::createPlane3(plane_trans);
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::offsetPoint(Point3SPtr point, Vector3SPtr dir, double offset) {
    Point3SPtr result;
#ifdef USE_CGAL
    Vector3 dir_normalized = *dir / CGAL::sqrt(dir->squared_length());
    Point3 p_moved = *point + (dir_normalized * offset);
#else
    Point3 p_moved = *point + (dir->normalize() * offset);
#endif
    result = KernelFactory::createPoint3(p_moved);
    DEBUG_SPTR(result);
    return result;
}

Vector3SPtr KernelWrapper::rotateVector(Vector3SPtr vector, Vector3SPtr axis, double angle) {
    Vector3SPtr result;
    Vector3SPtr v_n = KernelWrapper::normalize(axis);
    double n[3];
    for (unsigned int i = 0; i < 3; i++) {
        n[i] = (*v_n)[i];
    }
    double rotation[3][3];   // http://de.wikipedia.org/wiki/Drehmatrix
    rotation[0][0] = n[0]*n[0] * (1.0-cos(angle)) + cos(angle);
    rotation[0][1] = n[0]*n[1] * (1.0-cos(angle)) - n[2] * sin(angle);
    rotation[0][2] = n[0]*n[2] * (1.0-cos(angle)) + n[1] * sin(angle);
    rotation[1][0] = n[1]*n[0] * (1.0-cos(angle)) + n[2] * sin(angle);
    rotation[1][1] = n[1]*n[1] * (1.0-cos(angle)) + cos(angle);
    rotation[1][2] = n[1]*n[2] * (1.0-cos(angle)) - n[0] * sin(angle);
    rotation[2][0] = n[2]*n[0] * (1.0-cos(angle)) - n[1] * sin(angle);
    rotation[2][1] = n[2]*n[1] * (1.0-cos(angle)) + n[0] * sin(angle);
    rotation[2][2] = n[2]*n[2] * (1.0-cos(angle)) + cos(angle);
    double rotated[3];
    for (unsigned int r = 0; r < 3; r++) {
        rotated[r] = 0.0;
        for (unsigned int c = 0; c < 3; c++) {
            rotated[r] += rotation[r][c] * (*vector)[c];
        }
    }
    result = KernelFactory::createVector3(rotated[0], rotated[1], rotated[2]);
    DEBUG_SPTR(result);
    return result;
}

Plane3SPtr KernelWrapper::rotatePlane(Plane3SPtr plane, Line3SPtr line, double angle) {
    Plane3SPtr result;
    Point3SPtr point;
#ifdef USE_CGAL
    point = KernelFactory::createPoint3(line->point(0));
#else
    point = KernelFactory::createPoint3(line->point());
#endif
    Vector3SPtr dir = KernelFactory::createVector3(line);
    Vector3SPtr normal = KernelFactory::createVector3(plane);
    Vector3SPtr normal_rotated = rotateVector(normal, dir, angle);
    result = KernelFactory::createPlane3(point, normal_rotated);
    DEBUG_SPTR(result);
    return result;
}

int KernelWrapper::side(Plane3SPtr plane, Point3SPtr point) {
    int result = 0;
#ifdef USE_CGAL
    CGAL::Oriented_side side = plane->oriented_side(*point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
#else
    result = plane->side(*point);
#endif
    return result;
}

int KernelWrapper::orientation(Line3SPtr line1, Line3SPtr line2) {
    int result = 0;
#ifdef USE_CGAL
    Vector3 dir1 = line1->to_vector();
    Vector3 dir2 = line2->to_vector();
#else
    Vector3 dir1 = line1->direction();
    Vector3 dir2 = line2->direction();
#endif
    Point3 p0 = line1->point();
    Point3 p1 = p0 + dir1;
    Point3 p2 = line2->point();
    Plane3 plane(p0, p1, p2);
    Point3 point = p2 + dir2;
#ifdef USE_CGAL
    CGAL::Oriented_side side = plane.oriented_side(point);
    if (side == CGAL::ON_POSITIVE_SIDE) result = 1;
    if (side == CGAL::ON_NEGATIVE_SIDE) result = -1;
#else
    result = plane.side(point);
#endif
    return result;
}

double KernelWrapper::angle(Vector3SPtr v1, Vector3SPtr v2) {
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

double KernelWrapper::angle(Line3SPtr line1, Line3SPtr line2) {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(line1);
    Vector3SPtr v2 = KernelFactory::createVector3(line2);
    result = angle(v1, v2);
    return result;
}

double KernelWrapper::angle(Plane3SPtr plane, Line3SPtr line) {
    double result = 0.0;
    Vector3SPtr v_plane = KernelFactory::createVector3(plane);
    Vector3SPtr v_line = KernelFactory::createVector3(line);
    result = angle(v_plane, v_line);
    if (result > M_PI/2.0) {
        result = result - M_PI/2.0;
    } else {
        result = M_PI/2.0 - result;
    }
    return result;
}

double KernelWrapper::angle(Plane3SPtr plane1, Plane3SPtr plane2) {
    double result = 0.0;
    Vector3SPtr v1 = KernelFactory::createVector3(plane1);
    Vector3SPtr v2 = KernelFactory::createVector3(plane2);
    result = angle(v1, v2);
    return result;
}

bool KernelWrapper::isInside(Point3SPtr p, Point3SPtr p_box_1, Point3SPtr p_box_2) {
    bool result = true;
    for (unsigned int i = 0; i < 3; i++) {
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

Vector3SPtr KernelWrapper::cross(Vector3SPtr v1, Vector3SPtr v2) {
    Vector3SPtr result = Vector3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createVector3(CGAL::cross_product(*v1, *v2));
#else
    result = KernelFactory::createVector3(v1->cross(*v2));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::projection(Line3SPtr line, Point3SPtr point) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPoint3(line->projection(*point));
#else
    result = Point3SPtr(kernel::projection(&(*line), &(*point)));
#endif
    DEBUG_SPTR(result);
    return result;
}

Point3SPtr KernelWrapper::projection(Plane3SPtr plane, Point3SPtr point) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    result = KernelFactory::createPoint3(plane->projection(*point));
#else
    result = Point3SPtr(kernel::projection(&(*plane), &(*point)));
#endif
    DEBUG_SPTR(result);
    return result;
}

int KernelWrapper::comparePoints(Vector3SPtr v_dir, Point3SPtr p_1, Point3SPtr p_2) {
    int result = 0;
    double value = *v_dir * (*p_2 - *p_1);
    if (value > 0.0) {         // angle < M_PI/2.0
        result = -1;
    } else if (value < 0.0) {  // angle > M_PI/2.0
        result = 1;
    }
    return result;
}

Point3SPtr KernelWrapper::replaceCoord(Point3SPtr point, Point3SPtr replacement,
        unsigned int coord) {
    Point3SPtr result = Point3SPtr();
#ifdef USE_CGAL
    if (coord == 0) {
        result = Point3SPtr(new Point3(replacement->x(), point->y(), point->z()));
    } else if (coord == 1) {
        result = Point3SPtr(new Point3(point->x(), replacement->y(), point->z()));
    } else if (coord == 2) {
        result = Point3SPtr(new Point3(point->x(), point->y(), replacement->z()));
    }
#else
    if (coord == 0) {
        result = Point3SPtr(new Point3(replacement->getX(), point->getY(), point->getZ()));
    } else if (coord == 1) {
        result = Point3SPtr(new Point3(point->getX(), replacement->getY(), point->getZ()));
    } else if (coord == 2) {
        result = Point3SPtr(new Point3(point->getX(), point->getY(), replacement->getZ()));
    }
#endif
    DEBUG_SPTR(result);
    return result;
}

} }
