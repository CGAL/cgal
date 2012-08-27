// Copyright (c) 2001,2004  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion

#ifndef CGAL_FILTERED_KERNEL_H
#define CGAL_FILTERED_KERNEL_H

#include <CGAL/Filtered_kernel_fwd.h>
#include <CGAL/basic.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>

#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/internal/Exact_type_selector.h>

#include <CGAL/internal/Static_filters/Static_filters.h>

// This file contains the definition of a generic kernel filter.
//
// TODO:
// - at the moment, only the predicates are filtered.
//   Constructions will come later.
// - the kernel provides the traits interface, as well as type equality.
//   Having the global functions working is another story...
// - The converters are more a property of the types rather than anything else,
//   so maybe they should not be passed as template parameter, but use a
//   traits-like mechanism ?

namespace CGAL {

// CK = eventually rebound construction kernel (gets Point_2 from).
// Exact_kernel = exact kernel called when needed by the filter.
// Approximate_kernel = filtering "interval" kernel
template < typename CK >
struct Filtered_kernel_base
  : public CK
{
    typedef typename internal::Exact_type_selector<typename CK::RT>::Type  Exact_nt;
    typedef Simple_cartesian<Exact_nt>                           Exact_kernel;
    typedef Simple_cartesian<Interval_nt_advanced>               Approximate_kernel;
    typedef Cartesian_converter<CK, Exact_kernel>                C2E;
    typedef Cartesian_converter<CK, Approximate_kernel>          C2F;

    enum { Has_filtered_predicates = true };
    typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

    template < typename Kernel2 >
    struct Base {
        typedef typename CK::template Base<Kernel2> CK2;
        typedef Filtered_kernel_base<CK2>  Type;
    };

    template < typename T >
    struct Ambient_dimension {
        typedef typename T::Ambient_dimension type; // maybe not the right way...
    };

    template < typename T >
    struct Feature_dimension {
        typedef typename T::Feature_dimension type; // maybe not the right way...
    };

    // We change the predicates.
#define CGAL_Kernel_pred(P, Pf) \
    typedef Filtered_predicate<typename Exact_kernel::P, typename Approximate_kernel::P, C2E, C2F> P; \
    P Pf() const { return P(); }

    // We don't touch the constructions.
#define CGAL_Kernel_cons(Y,Z)

#include <CGAL/Kernel/interface_macros.h>

};

template < typename CK >
struct Static_filters_base
  : public internal::Static_filters< Filtered_kernel_base<CK> >
{
    template < typename Kernel2 >
    struct Base {
        typedef typename CK::template Base<Kernel2>::Type  CK2;
        typedef Static_filters_base<CK2>                   Type;
    };
};

#ifdef CGAL_NO_STATIC_FILTERS
template < typename CK, bool UseStaticFilters = false >
#else
template < typename CK, bool UseStaticFilters = true >
#endif
struct Filtered_kernel_adaptor
  : public Filtered_kernel_base<CK>
{
	enum { Has_static_filters = false };
};

template < typename CK >
struct Filtered_kernel_adaptor<CK, true>
  : public Static_filters_base<CK>
{
	enum { Has_static_filters = true };
};

// UseStaticFilters has a default value, depending on
// CGAL_NO_STATIC_FILTERS. See in <CGAL/Filtered_kernel_fwd.h>.
template < typename CK, bool UseStaticFilters >
struct Filtered_kernel
  : public Filtered_kernel_adaptor<
               Type_equality_wrapper<
                   typename CK:: template Base< Filtered_kernel<CK, UseStaticFilters> >::Type,
                   Filtered_kernel<CK, UseStaticFilters> >,
	       UseStaticFilters >
{};

} //namespace CGAL

#endif // CGAL_FILTERED_KERNEL_H
