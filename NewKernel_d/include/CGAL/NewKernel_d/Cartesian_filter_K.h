// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_CARTESIAN_FILTER_K_H
#define CGAL_KERNEL_D_CARTESIAN_FILTER_K_H

#include <CGAL/basic.h>
#include <CGAL/NewKernel_d/KernelD_converter.h>
#include <CGAL/NewKernel_d/Filtered_predicate2.h>
#include <boost/mpl/if.hpp>
#include <boost/mpl/and.hpp>

namespace CGAL {

template < typename Base_, typename AK_, typename EK_ >
struct Cartesian_filter_K : public Base_,
  private Store_kernel<AK_>, private Store_kernel2<EK_>
{
    constexpr Cartesian_filter_K(){}
    constexpr Cartesian_filter_K(int d):Base_(d){}
    //FIXME: or do we want an instance of AK and EK belonging to this kernel,
    //instead of a reference to external ones?
    constexpr Cartesian_filter_K(AK_ const&a,EK_ const&b):Base_(),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    constexpr Cartesian_filter_K(int d,AK_ const&a,EK_ const&b):Base_(d),Store_kernel<AK_>(a),Store_kernel2<EK_>(b){}
    typedef Base_ Kernel_base;
    typedef AK_ AK;
    typedef EK_ EK;
    typedef typename Store_kernel<AK_>::reference_type AK_rt;
    AK_rt approximate_kernel()const{return this->kernel();}
    typedef typename Store_kernel2<EK_>::reference2_type EK_rt;
    EK_rt exact_kernel()const{return this->kernel2();}

    // MSVC is too dumb to perform the empty base optimization.
    typedef boost::mpl::and_<
      internal::Do_not_store_kernel<Kernel_base>,
      internal::Do_not_store_kernel<AK>,
      internal::Do_not_store_kernel<EK> > Do_not_store_kernel;

    //TODO: C2A/C2E could be able to convert *this into this->kernel() or this->kernel2().
    typedef KernelD_converter<Kernel_base,AK> C2A;
    typedef KernelD_converter<Kernel_base,EK> C2E;

    // fix the types
    // TODO: only fix some types, based on some criterion?
    template<class T> struct Type : Get_type<Kernel_base,T> {};

    template<class T,class D=void,class=typename Get_functor_category<Cartesian_filter_K,T>::type> struct Functor :
      Inherit_functor<Kernel_base,T,D> {};
    template<class T,class D> struct Functor<T,D,Predicate_tag> {
      typedef typename Get_functor<AK, T>::type AP;
      typedef typename Get_functor<EK, T>::type EP;
      typedef Filtered_predicate2<EP,AP,C2E,C2A> type;
    };
// TODO:
//    template<class T> struct Functor<T,No_filter_tag,Predicate_tag> :
//	    Kernel_base::template Functor<T,No_filter_tag> {};
// TODO:
// detect when Less_cartesian_coordinate doesn't need filtering
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_FILTER_K_H
