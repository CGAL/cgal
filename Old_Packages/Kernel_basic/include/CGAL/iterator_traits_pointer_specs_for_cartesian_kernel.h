// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 
//
// file          : iterator_traits_pointer_specs_for_cartesian_kernel.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra <Stefan.Schirra@mpi-sb.mpg.de>
//
// coordinator   : MPI, Saarbruecken
// ======================================================================


#ifndef CGAL_ITERATOR_TRAITS_POINTER_SPECS_FOR_CARTESIAN_KERNEL_H
#define CGAL_ITERATOR_TRAITS_POINTER_SPECS_FOR_CARTESIAN_KERNEL_H

// to be included in Cartesian.h

#ifdef CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#include <CGAL/user_classes.h>
#include <CGAL/cartesian_classes.h>

#define CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C(NT)                                              \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::Point_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::Point_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Point_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::Point_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::Point_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Point_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::Vector_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Vector_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::Vector_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Vector_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::Vector_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Vector_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::Vector_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Vector_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                 \
    template <>                                                                       \
    struct iterator_traits<const CGAL::Direction_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::Direction_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef const CGAL::Direction_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Direction_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                \
    template <>                                                                       \
    struct iterator_traits<CGAL::Direction_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::Direction_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef CGAL::Direction_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Direction_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                \
__STL_END_NAMESPACE                                                                   \
__STL_BEGIN_NAMESPACE                                                            \
    template <>                                                                  \
    struct iterator_traits<const CGAL::Line_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::Line_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef const CGAL::Line_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Line_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                           \
    template <>                                                                  \
    struct iterator_traits<CGAL::Line_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::Line_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef CGAL::Line_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Line_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                           \
__STL_END_NAMESPACE                                                              \
__STL_BEGIN_NAMESPACE                                                               \
    template <>                                                                     \
    struct iterator_traits<const CGAL::Segment_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::Segment_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef const CGAL::Segment_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Segment_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                              \
    template <>                                                                     \
    struct iterator_traits<CGAL::Segment_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::Segment_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef CGAL::Segment_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Segment_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                              \
__STL_END_NAMESPACE                                                                 \
__STL_BEGIN_NAMESPACE                                                           \
    template <>                                                                 \
    struct iterator_traits<const CGAL::Ray_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::Ray_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef const CGAL::Ray_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Ray_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                          \
    template <>                                                                 \
    struct iterator_traits<CGAL::Ray_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::Ray_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef CGAL::Ray_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Ray_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                          \
__STL_END_NAMESPACE                                                             \
__STL_BEGIN_NAMESPACE                                                                     \
    template <>                                                                           \
    struct iterator_traits<const CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                             iterator_category; \
        typedef CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                              difference_type;   \
        typedef const CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                    \
    template <>                                                                           \
    struct iterator_traits<CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                             iterator_category; \
        typedef CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                              difference_type;   \
        typedef CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Iso_rectangle_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                    \
__STL_END_NAMESPACE                                                                       \
__STL_BEGIN_NAMESPACE                                                                \
    template <>                                                                      \
    struct iterator_traits<const CGAL::Triangle_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::Triangle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef const CGAL::Triangle_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Triangle_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                               \
    template <>                                                                      \
    struct iterator_traits<CGAL::Triangle_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::Triangle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef CGAL::Triangle_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Triangle_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                               \
__STL_END_NAMESPACE                                                                  \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::Circle_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Circle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::Circle_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Circle_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::Circle_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Circle_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::Circle_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Circle_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                          \
    template <>                                                                                \
    struct iterator_traits<const CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef const CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                         \
    template <>                                                                                \
    struct iterator_traits<CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Aff_transformation_2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                         \
__STL_END_NAMESPACE                                                                            \

#define CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C(NT)                                              \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::Point_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::Point_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Point_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::Point_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::Point_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Point_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::Vector_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Vector_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::Vector_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Vector_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::Vector_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::Vector_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::Vector_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Vector_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                 \
    template <>                                                                       \
    struct iterator_traits<const CGAL::Direction_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::Direction_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef const CGAL::Direction_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Direction_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                \
    template <>                                                                       \
    struct iterator_traits<CGAL::Direction_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::Direction_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef CGAL::Direction_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Direction_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                \
__STL_END_NAMESPACE                                                                   \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::Plane_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Plane_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::Plane_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Plane_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::Plane_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Plane_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::Plane_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Plane_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                            \
    template <>                                                                  \
    struct iterator_traits<const CGAL::Line_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::Line_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef const CGAL::Line_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Line_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                           \
    template <>                                                                  \
    struct iterator_traits<CGAL::Line_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::Line_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef CGAL::Line_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Line_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                           \
__STL_END_NAMESPACE                                                              \
__STL_BEGIN_NAMESPACE                                                               \
    template <>                                                                     \
    struct iterator_traits<const CGAL::Segment_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::Segment_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef const CGAL::Segment_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Segment_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                              \
    template <>                                                                     \
    struct iterator_traits<CGAL::Segment_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::Segment_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef CGAL::Segment_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Segment_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                              \
__STL_END_NAMESPACE                                                                 \
__STL_BEGIN_NAMESPACE                                                           \
    template <>                                                                 \
    struct iterator_traits<const CGAL::Ray_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::Ray_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef const CGAL::Ray_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Ray_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                          \
    template <>                                                                 \
    struct iterator_traits<CGAL::Ray_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::Ray_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef CGAL::Ray_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Ray_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                          \
__STL_END_NAMESPACE                                                             \
__STL_BEGIN_NAMESPACE                                                                \
    template <>                                                                      \
    struct iterator_traits<const CGAL::Triangle_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::Triangle_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef const CGAL::Triangle_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Triangle_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                               \
    template <>                                                                      \
    struct iterator_traits<CGAL::Triangle_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::Triangle_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef CGAL::Triangle_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Triangle_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                               \
__STL_END_NAMESPACE                                                                  \
__STL_BEGIN_NAMESPACE                                                                   \
    template <>                                                                         \
    struct iterator_traits<const CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                           iterator_category; \
        typedef CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                            difference_type;   \
        typedef const CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                  \
    template <>                                                                         \
    struct iterator_traits<CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                           iterator_category; \
        typedef CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                            difference_type;   \
        typedef CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Tetrahedron_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                  \
__STL_END_NAMESPACE                                                                     \
__STL_BEGIN_NAMESPACE                                                                          \
    template <>                                                                                \
    struct iterator_traits<const CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef const CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                         \
    template <>                                                                                \
    struct iterator_traits<CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Aff_transformation_3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                         \
__STL_END_NAMESPACE                                                                            \

#define CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC(NT)                                              \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::Point_d< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_d< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::Point_d< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Point_d< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::Point_d< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::Point_d< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::Point_d< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Point_d< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \

#define CGAL_ITERATOR_TRAITS_POINTER_SPECC2(NT)                                               \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::PointC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::PointC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::PointC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::PointC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::PointC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::PointC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::VectorC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::VectorC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::VectorC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::VectorC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::VectorC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::VectorC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::VectorC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::VectorC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                 \
    template <>                                                                       \
    struct iterator_traits<const CGAL::DirectionC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::DirectionC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef const CGAL::DirectionC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::DirectionC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                \
    template <>                                                                       \
    struct iterator_traits<CGAL::DirectionC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::DirectionC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef CGAL::DirectionC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::DirectionC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                \
__STL_END_NAMESPACE                                                                   \
__STL_BEGIN_NAMESPACE                                                            \
    template <>                                                                  \
    struct iterator_traits<const CGAL::LineC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::LineC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef const CGAL::LineC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::LineC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                           \
    template <>                                                                  \
    struct iterator_traits<CGAL::LineC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::LineC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef CGAL::LineC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::LineC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                           \
__STL_END_NAMESPACE                                                              \
__STL_BEGIN_NAMESPACE                                                               \
    template <>                                                                     \
    struct iterator_traits<const CGAL::SegmentC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::SegmentC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef const CGAL::SegmentC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::SegmentC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                              \
    template <>                                                                     \
    struct iterator_traits<CGAL::SegmentC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::SegmentC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef CGAL::SegmentC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::SegmentC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                              \
__STL_END_NAMESPACE                                                                 \
__STL_BEGIN_NAMESPACE                                                           \
    template <>                                                                 \
    struct iterator_traits<const CGAL::RayC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::RayC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef const CGAL::RayC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::RayC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                          \
    template <>                                                                 \
    struct iterator_traits<CGAL::RayC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::RayC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef CGAL::RayC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::RayC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                          \
__STL_END_NAMESPACE                                                             \
__STL_BEGIN_NAMESPACE                                                                     \
    template <>                                                                           \
    struct iterator_traits<const CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                             iterator_category; \
        typedef CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                              difference_type;   \
        typedef const CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                    \
    template <>                                                                           \
    struct iterator_traits<CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                             iterator_category; \
        typedef CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                              difference_type;   \
        typedef CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Iso_rectangleC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                    \
__STL_END_NAMESPACE                                                                       \
__STL_BEGIN_NAMESPACE                                                                \
    template <>                                                                      \
    struct iterator_traits<const CGAL::TriangleC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::TriangleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef const CGAL::TriangleC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::TriangleC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                               \
    template <>                                                                      \
    struct iterator_traits<CGAL::TriangleC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::TriangleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef CGAL::TriangleC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::TriangleC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                               \
__STL_END_NAMESPACE                                                                  \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::CircleC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::CircleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::CircleC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::CircleC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::CircleC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::CircleC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::CircleC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::CircleC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                          \
    template <>                                                                                \
    struct iterator_traits<const CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef const CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                         \
    template <>                                                                                \
    struct iterator_traits<CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Aff_transformationC2< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                         \
__STL_END_NAMESPACE                                                                            \

#define CGAL_ITERATOR_TRAITS_POINTER_SPECC3(NT)                                               \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::PointC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::PointC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::PointC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::PointC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::PointC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::PointC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                              \
    template <>                                                                    \
    struct iterator_traits<const CGAL::VectorC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::VectorC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef const CGAL::VectorC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::VectorC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                             \
    template <>                                                                    \
    struct iterator_traits<CGAL::VectorC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                      iterator_category; \
        typedef CGAL::VectorC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                       difference_type;   \
        typedef CGAL::VectorC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::VectorC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                             \
__STL_END_NAMESPACE                                                                \
__STL_BEGIN_NAMESPACE                                                                 \
    template <>                                                                       \
    struct iterator_traits<const CGAL::DirectionC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::DirectionC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef const CGAL::DirectionC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::DirectionC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                \
    template <>                                                                       \
    struct iterator_traits<CGAL::DirectionC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                         iterator_category; \
        typedef CGAL::DirectionC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                          difference_type;   \
        typedef CGAL::DirectionC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::DirectionC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                \
__STL_END_NAMESPACE                                                                   \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::PlaneC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PlaneC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::PlaneC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::PlaneC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::PlaneC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PlaneC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::PlaneC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::PlaneC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
__STL_BEGIN_NAMESPACE                                                            \
    template <>                                                                  \
    struct iterator_traits<const CGAL::LineC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::LineC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef const CGAL::LineC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::LineC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                           \
    template <>                                                                  \
    struct iterator_traits<CGAL::LineC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                    iterator_category; \
        typedef CGAL::LineC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                     difference_type;   \
        typedef CGAL::LineC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::LineC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                           \
__STL_END_NAMESPACE                                                              \
__STL_BEGIN_NAMESPACE                                                               \
    template <>                                                                     \
    struct iterator_traits<const CGAL::SegmentC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::SegmentC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef const CGAL::SegmentC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::SegmentC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                              \
    template <>                                                                     \
    struct iterator_traits<CGAL::SegmentC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                       iterator_category; \
        typedef CGAL::SegmentC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                        difference_type;   \
        typedef CGAL::SegmentC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::SegmentC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                              \
__STL_END_NAMESPACE                                                                 \
__STL_BEGIN_NAMESPACE                                                           \
    template <>                                                                 \
    struct iterator_traits<const CGAL::RayC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::RayC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef const CGAL::RayC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::RayC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                          \
    template <>                                                                 \
    struct iterator_traits<CGAL::RayC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                   iterator_category; \
        typedef CGAL::RayC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                    difference_type;   \
        typedef CGAL::RayC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::RayC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                          \
__STL_END_NAMESPACE                                                             \
__STL_BEGIN_NAMESPACE                                                                \
    template <>                                                                      \
    struct iterator_traits<const CGAL::TriangleC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::TriangleC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef const CGAL::TriangleC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::TriangleC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                               \
    template <>                                                                      \
    struct iterator_traits<CGAL::TriangleC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                        iterator_category; \
        typedef CGAL::TriangleC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                         difference_type;   \
        typedef CGAL::TriangleC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::TriangleC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                               \
__STL_END_NAMESPACE                                                                  \
__STL_BEGIN_NAMESPACE                                                                   \
    template <>                                                                         \
    struct iterator_traits<const CGAL::TetrahedronC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                           iterator_category; \
        typedef CGAL::TetrahedronC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                            difference_type;   \
        typedef const CGAL::TetrahedronC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::TetrahedronC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                  \
    template <>                                                                         \
    struct iterator_traits<CGAL::TetrahedronC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                           iterator_category; \
        typedef CGAL::TetrahedronC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                            difference_type;   \
        typedef CGAL::TetrahedronC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::TetrahedronC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                  \
__STL_END_NAMESPACE                                                                     \
__STL_BEGIN_NAMESPACE                                                                          \
    template <>                                                                                \
    struct iterator_traits<const CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef const CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                                         \
    template <>                                                                                \
    struct iterator_traits<CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                                  iterator_category; \
        typedef CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                                   difference_type;   \
        typedef CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::Aff_transformationC3< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                                         \
__STL_END_NAMESPACE                                                                            \

#define CGAL_ITERATOR_TRAITS_POINTER_SPECCD(NT)                                   \
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::PointCd< NT >*> {                          \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointCd< NT >         value_type;                           \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::PointCd< NT >*  pointer;                              \
        typedef const CGAL::PointCd< NT >&  reference;                            \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::PointCd< NT >*> {                                \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointCd< NT >         value_type;                           \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::PointCd< NT >*        pointer;                              \
        typedef CGAL::PointCd< NT >&        reference;                            \
    };                                                                            \
__STL_END_NAMESPACE                                                               \

/* 
__STL_BEGIN_NAMESPACE                                                             \
    template <>                                                                   \
    struct iterator_traits<const CGAL::PointCd< CGAL::Cartesian< NT > >*> {       \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointCd< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef const CGAL::PointCd< CGAL::Cartesian< NT > >*  pointer;           \
        typedef const CGAL::PointCd< CGAL::Cartesian< NT > >&  reference;         \
    };                                                                            \
    template <>                                                                   \
    struct iterator_traits<CGAL::PointCd< CGAL::Cartesian< NT > >*> {             \
        typedef random_access_iterator_tag                     iterator_category; \
        typedef CGAL::PointCd< CGAL::Cartesian< NT > >         value_type;        \
        typedef ptrdiff_t                                      difference_type;   \
        typedef CGAL::PointCd< CGAL::Cartesian< NT > >*        pointer;           \
        typedef CGAL::PointCd< CGAL::Cartesian< NT > >&        reference;         \
    };                                                                            \
__STL_END_NAMESPACE                                                               \
*/



CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( int )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( int )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( int )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( long )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( long )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( long )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( float )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( float )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( float )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( double )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( double )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( double )

CGAL_ITERATOR_TRAITS_POINTER_SPECC2( int )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( int )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( int )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( long )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( long )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( long )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( float )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( float )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( float )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( double )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( double )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( double )

class leda_real;
class leda_integer;
class leda_rational;
class leda_bigfloat;

CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( leda_bigfloat )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( leda_bigfloat )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( leda_bigfloat )

CGAL_ITERATOR_TRAITS_POINTER_SPECC2( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( leda_real )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( leda_integer )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( leda_rational )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( leda_bigfloat )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( leda_bigfloat )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( leda_bigfloat )

namespace CGAL { class Gmpz; }

CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Gmpz )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Gmpz )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Gmpz )

CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Gmpz )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Gmpz )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Gmpz )


#ifdef CGAL_QUOTIENT_H
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_2C( CGAL::Quotient<CGAL::Gmpz> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_3C( CGAL::Quotient<CGAL::Gmpz> )
CGAL_ITERATOR_TRAITS_POINTER_SPEC_DC( CGAL::Quotient<CGAL::Gmpz> )

CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<int> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<long> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<float> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<double> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<leda_real> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<leda_integer> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC2( CGAL::Quotient<CGAL::Gmpz> )
CGAL_ITERATOR_TRAITS_POINTER_SPECC3( CGAL::Quotient<CGAL::Gmpz> )
CGAL_ITERATOR_TRAITS_POINTER_SPECCD( CGAL::Quotient<CGAL::Gmpz> )

#endif // CGAL_QUOTIENT_H

#endif // CGAL_LIMITED_ITERATOR_TRAITS_SUPPORT
#endif // CGAL_ITERATOR_TRAITS_POINTER_SPECS_FOR_CARTESIAN_KERNEL_H
