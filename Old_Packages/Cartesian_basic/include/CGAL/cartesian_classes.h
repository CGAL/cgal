// revision      : $Revision$
// revision_date : $Date$
// authors       : Herve Bronnimann

#ifndef CGAL_CARTESIAN_CLASSES_H
#define CGAL_CARTESIAN_CLASSES_H

#include <CGAL/basic_classes.h>

#ifndef CGAL_CFG_NO_ADVANCED_KERNEL

CGAL_BEGIN_NAMESPACE

// The following scheme is proposed by Michael Hoffmann
// It uses partial specialization and default template parameters
// and a curiously recurring pattern
// It allows for a fully extendible kernel

template < class R, class T = typename R::Rep_tag > class Point_2;
template < class R, class T = typename R::Rep_tag > class Vector_2;
template < class R, class T = typename R::Rep_tag > class Direction_2;
template < class R, class T = typename R::Rep_tag > class Line_2;
template < class R, class T = typename R::Rep_tag > class Ray_2;
template < class R, class T = typename R::Rep_tag > class Segment_2;
template < class R, class T = typename R::Rep_tag > class Triangle_2;
template < class R, class T = typename R::Rep_tag > class Circle_2;
template < class R, class T = typename R::Rep_tag > class Data_accessor_2;
template < class PT, class DA > class ConicCPA2;
template < class R, class T = typename R::Rep_tag > class Iso_rectangle_2;
template < class R, class T = typename R::Rep_tag > class Aff_transformation_2;

template < class R, class T = typename R::Rep_tag > class Plane_3;
template < class R, class T = typename R::Rep_tag > class Point_3;
template < class R, class T = typename R::Rep_tag > class Vector_3;
template < class R, class T = typename R::Rep_tag > class Direction_3;
template < class R, class T = typename R::Rep_tag > class Line_3;
template < class R, class T = typename R::Rep_tag > class Ray_3;
template < class R, class T = typename R::Rep_tag > class Segment_3;
template < class R, class T = typename R::Rep_tag > class Triangle_3;
template < class R, class T = typename R::Rep_tag > class Tetrahedron_3;
template < class R, class T = typename R::Rep_tag > class Iso_cuboid_3;
template < class R, class T = typename R::Rep_tag > class Aff_transformation_3;

template < class R > class PointCd;

CGAL_END_NAMESPACE

#else 

CGAL_BEGIN_NAMESPACE

// The following scheme is proposed by Stefan Schirra
// It does not use partial specialization
// It is a partially extendible kernel, but you have to redefine and
// reinstantiate all the member classes in a new kernel
// even if only one differs

template < class R > class PointC2;
template < class R > class VectorC2;
template < class R > class DirectionC2;
template < class R > class LineC2;
template < class R > class RayC2;
template < class R > class SegmentC2;
template < class R > class TriangleC2;
template < class R > class CircleC2;
template < class R > class Data_accessorC2;
template < class PT, class DA > class ConicCPA2;
template < class R > class Iso_rectangleC2;
template < class R > class Aff_transformationC2;

template < class R > class PlaneC3;
template < class R > class PointC3;
template < class R > class VectorC3;
template < class R > class DirectionC3;
template < class R > class LineC3;
template < class R > class RayC3;
template < class R > class SegmentC3;
template < class R > class TriangleC3;
template < class R > class TetrahedronC3;
template < class R > class Iso_cuboidC3;
template < class R > class Aff_transformationC3;

template < class R > class PointCd;

CGAL_END_NAMESPACE

// We also need the wrapper classes Point_2<R> etc.
// We include them (they are common to Cartesian and Homogeneous)
#include <CGAL/user_classes.h>

#endif // CGAL_CFG_NO_ADVANCED_KERNEL

#endif // CGAL_CARTESIAN_CLASSES_H
