// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Filtered_kernel.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis
//
// ============================================================================

#ifndef CGAL_FILTERED_KERNEL_H
#define CGAL_FILTERED_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

// This file contains the definition of a generic kernel filter.
//
// TODO:
// - at the moment, it's restricted to IA filtering, but this should be
//   generalized to allow other kinds of filters (static...).
// - at the moment, only the predicates are filtered.
//   Constructions will come later.
// - the kernel provides the traits interface, as well as type equality.
//   Having the global functions working is another story...
// - The converters are more a property of the types rather than anything else,
//   so maybe they should not be passed as template parameter, but use a
//   traits-like mecanism ?

CGAL_BEGIN_NAMESPACE

// CK = construction kernel.
// EK = exact kernel called when needed by the filter.
// FK = filtering kernel
template < typename CK, typename Kernel >
class Filtered_kernel_base
  : public CK::template base<Kernel>::other
{
    typedef typename CK::template base<Kernel>::other   Base;
    // Hardcoded for now.
    typedef Simple_cartesian<Quotient<MP_Float> >               EK;
    typedef Simple_cartesian<Interval_nt_advanced>   FK;
    typedef Cartesian_converter<Base, EK>            C2E;
    typedef Cartesian_converter<Base, FK,
                                Interval_converter<typename Base::RT> > C2F;
public:

    template < typename Kernel2 >
    struct base { typedef Filtered_kernel_base<CK, Kernel2>  other; };

    // What to do with the tag ?
    // Probably this should not exist, should it ?
    // struct filter_tag{};
    // typedef filter_tag                                     Kernel_tag;
    // typedef typename CK::Kernel_tag                       Kernel_tag;
    // typedef typename CK::Rep_tag                          Rep_tag;

    // We change the predicates.
#define CGAL_Kernel_pred(P, Pf) \
    typedef Filtered_predicate<typename EK::P, typename FK::P, \
	                     C2E, C2F> P; \
    P Pf() const { return P(); }

    // We don't touch the constructions.
#define CGAL_Kernel_cons(Y,Z)

#include <CGAL/Kernel/interface_macros.h>

};

template <class CK>
struct Filtered_kernel
  : public Type_equality_wrapper<
                Filtered_kernel_base< CK, Filtered_kernel<CK> >,
                Filtered_kernel<CK> >
{
  template < typename Kernel >
  struct base { typedef Filtered_kernel_base<CK, Kernel>   other; };
};

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_KERNEL_H
