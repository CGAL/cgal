// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Simple_homogeneous.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra, Sylvain Pion
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_SIMPLE_HOMOGENEOUS_H
#define CGAL_SIMPLE_HOMOGENEOUS_H

#include <CGAL/Homogeneous/Homogeneous_base.h>
#include <CGAL/Simple_Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Kernel/function_objects.h>

#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

template < typename RT_, typename FT_ = Quotient<RT_> >
class Simple_homogeneous
  : public Type_equality_wrapper<
              Homogeneous_base< Simple_homogeneous<RT_, FT_> >,
              Simple_homogeneous<RT_, FT_> >
{
    typedef Simple_homogeneous<RT_, FT_>                  Kernel;
public:

    typedef RT_                                           RT;
    typedef FT_                                           FT;

    // The mecanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle {
        typedef Simple_Handle_for<T>    type;
    };

    // Undocumented stuff.
    typedef Data_accessorH2<Kernel>                       Data_accessor_2;
    typedef CGAL::Conic_2<Kernel>                         Conic_2;

    // TODO: cleanup (use Rational_traits<> instead)
    static FT make_FT(const RT & num, const RT& denom)
    { return FT(num, denom); }

    static FT make_FT(const RT & num)
    { return FT(num); }

    static RT FT_numerator(const FT &r)
    { return r.numerator(); }

    static RT FT_denominator(const FT &r)
    { return r.denominator(); }

    // Functors types and access functions.
#define CGAL_Kernel_pred(Y,Z) typedef CGALi::Y<Kernel> Y; \
                              Y Z() const { return Y(); }
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

#include <CGAL/Kernel/interface_macros.h>

};

CGAL_END_NAMESPACE

CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(CGAL::Simple_homogeneous)

#endif // CGAL_SIMPLE_HOMOGENEOUS_H
