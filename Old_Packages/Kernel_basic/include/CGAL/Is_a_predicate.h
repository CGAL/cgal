// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Is_a_predicate.h
// package       : Kernel_basic
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================

#ifndef CGAL_IS_A_PREDICATE_H
#define CGAL_IS_A_PREDICATE_H

// How to determine if a kernel functor is a predicate or a construction.

#include <CGAL/basic.h>
#include <CGAL/enum.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

// By default it's a construction
template <typename Return_type>
struct Return_type_of_predicate {
    typedef CGAL::Tag_false type;
};

// Specializations for predicates
template <>
struct Return_type_of_predicate<bool> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Sign> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Orientation> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Oriented_side> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Bounded_side> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Comparison_result> {
    typedef CGAL::Tag_true type;
};

template <>
struct Return_type_of_predicate<CGAL::Angle> {
    typedef CGAL::Tag_true type;
};

} // namespace CGALi

template <typename Functor>
struct Is_a_predicate {
  typedef typename CGALi::Return_type_of_predicate<
                   typename Functor::result_type>::type type;
};

CGAL_END_NAMESPACE

#endif // CGAL_IS_A_PREDICATE_H
