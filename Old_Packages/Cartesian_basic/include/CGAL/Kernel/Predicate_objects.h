// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_KERNEL_PREDICATE_OBJECTS_H
#define CGAL_KERNEL_PREDICATE_OBJECTS_H

#include <CGAL/Kernel/Predicate_objects_2.h>
#include <CGAL/Kernel/Predicate_objects_3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Kernel_predicate_objects
{
public:
typedef typename R::FT                         FT;
typedef typename R::RT                         RT;
typedef typename R::Point_3                    Point_3;
typedef typename R::Vector_3                   Vector_3;
typedef typename R::Direction_3                Direction_3;
typedef typename R::Line_3                     Line_3;
typedef typename R::Plane_3                    Plane_3;
typedef typename R::Ray_3                      Ray_3;
typedef typename R::Segment_3                  Segment_3;
typedef typename R::Triangle_3                 Triangle_3;
typedef typename R::Tetrahedron_3              Tetrahedron_3;
typedef typename R::Aff_transformation_3       Aff_transformation_3;

CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_2(typename Kernel_predicate_objects_2<R>)
CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_3(typename Kernel_predicate_objects_3<R>) 
};

CGAL_END_NAMESPACE

// This macro is provided for convenience in defining the Kernel
// function objects inside a new representation class.
// See Cartesian_3.h and Cartesian.h

#define CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS(PR) \
CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_2(PR) \
CGAL_UNPACK_KERNEL_PREDICATE_OBJECTS_3(PR) 

#endif // CGAL_KERNEL_PREDICATE_OBJECTS_3_H
