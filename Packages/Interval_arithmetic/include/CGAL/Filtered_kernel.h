// Copyright (c) 2001,2004  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_FILTERED_KERNEL_H
#define CGAL_FILTERED_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>

#ifndef CGAL_NO_STATIC_FILTERS
#  include <CGAL/Static_filters.h>
#endif

// This file contains the definition of a generic kernel filter.
//
// TODO:
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
  : public CK::template Base<Kernel>::Type
{
    typedef typename CK::template Base<Kernel>::Type   Kernel_base;
public:

    // Hardcoded for now.
    typedef Simple_cartesian<Quotient<MP_Float> >    EK;
    typedef Simple_cartesian<Interval_nt_advanced>   FK;
    typedef Cartesian_converter<Kernel, EK,
                                // we need to specify the default arg otherwise
                                // Kernel would be instantiated.
                                NT_converter<typename Kernel_base::RT,
                                             typename EK::RT> >     C2E;
    typedef Cartesian_converter<Kernel, FK,
                                To_interval<typename Kernel_base::RT> > C2F;

    template < typename Kernel2 >
    struct Base { typedef Filtered_kernel_base<CK, Kernel2>  Type; };

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

#ifndef CGAL_NO_STATIC_FILTERS
template < typename CK, typename Kernel >
class Static_filters_base
  : public Static_filters<Filtered_kernel_base<CK, Kernel> >
{
    template < typename Kernel2 >
    struct Base { typedef Static_filters_base<CK, Kernel2>  Type; };
};
#endif

template <class CK>
struct Filtered_kernel_adaptor
#ifndef CGAL_NO_STATIC_FILTERS
  : public Static_filters_base< CK, Filtered_kernel_adaptor<CK> >
#else
  : public Filtered_kernel_base< CK, Filtered_kernel_adaptor<CK> >
#endif
{};

template <class CK>
struct Filtered_kernel_without_type_equality
#ifndef CGAL_NO_STATIC_FILTERS
  : public Static_filters_base< CK, Filtered_kernel_without_type_equality<CK> >
#else
  : public Filtered_kernel_base< CK, Filtered_kernel_without_type_equality<CK> >
#endif
{};

template <class CK>
struct Filtered_kernel
  : public Type_equality_wrapper< 
#ifndef CGAL_NO_STATIC_FILTERS
             Static_filters_base< CK, Filtered_kernel<CK> >,
#else
             Filtered_kernel_base< CK, Filtered_kernel<CK> >,
#endif
             Filtered_kernel<CK> >
{};

CGAL_END_NAMESPACE

#endif // CGAL_FILTERED_KERNEL_H
