// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann

#ifndef CGAL_KERNEL_CONSTRUCTION_OBJECTS_H
#define CGAL_KERNEL_CONSTRUCTION_OBJECTS_H

#include <CGAL/Kernel/Construction_objects_2.h>
#include <CGAL/Kernel/Construction_objects_3.h>

CGAL_BEGIN_NAMESPACE

template < class R >
class Kernel_construction_objects
{
public:
typedef typename R::FT                      FT;
typedef typename R::RT                      RT;
typedef typename R::Point_2                 Point_2;
typedef typename R::Vector_2                Vector_2;
typedef typename R::Direction_2             Direction_2;
typedef typename R::Segment_2               Segment_2;
typedef typename R::Line_2                  Line_2;
typedef typename R::Ray_2                   Ray_2;
typedef typename R::Triangle_2              Triangle_2;
typedef typename R::Circle_2                Circle_2;
typedef typename R::Iso_rectangle_2         Iso_rectangle_2;
typedef typename R::Aff_transformation_2    Aff_transformation_2;

CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS_2(typename 
                                      Kernel_construction_objects_2<R>)
CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS_3(typename 
                                      Kernel_construction_objects_3<R>)
};

CGAL_END_NAMESPACE

// This macro is provided for convenience in defining the Kernel
// function objects inside a new representation class.
// See Cartesian_2.h and Cartesian.h

#define CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS(CO) \
CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS_2(CO) \
CGAL_UNPACK_KERNEL_CONSTRUCTION_OBJECTS_3(CO)

#endif // CGAL_KERNEL_CONSTRUCTION_OBJECTS_2_H
