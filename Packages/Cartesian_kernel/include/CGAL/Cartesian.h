// ======================================================================
//
// Copyright (c) 2000,2001,2002,2003 The CGAL Consortium
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
// file          : include/CGAL/Cartesian.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================

#ifndef CGAL_CARTESIAN_H
#define CGAL_CARTESIAN_H

#include <CGAL/Cartesian/Cartesian_base.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>

CGAL_BEGIN_NAMESPACE

template < typename FT_, typename Kernel >
struct Cartesian_base_ref_count
  : public Cartesian_base< Kernel >
{
    typedef FT_                                           RT;
    typedef FT_                                           FT;

    // The mecanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef Handle_for<T>    type; };

    template < typename Kernel2 >
    struct base { typedef Cartesian_base_ref_count<FT_, Kernel2>  other; };

    // TODO: cleanup (use Rational_traits<> instead)
    static   FT make_FT(const RT & num, const RT& denom) { return num/denom;}
    static   FT make_FT(const RT & num)                  { return num;}
    static   RT FT_numerator(const FT &r)                { return r;}
    static   RT FT_denominator(const FT &)               { return RT(1);}
};

template < typename FT_ >
struct Cartesian
  : public Type_equality_wrapper<
                Cartesian_base_ref_count<FT_, Cartesian<FT_> >,
                Cartesian<FT_> >
{};

CGAL_END_NAMESPACE

CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(CGAL::Cartesian)

#endif // CGAL_CARTESIAN_H
