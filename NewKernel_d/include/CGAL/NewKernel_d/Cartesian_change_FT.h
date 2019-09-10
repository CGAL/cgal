// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Marc Glisse

#ifndef CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H
#define CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/NewKernel_d/Cartesian_complete.h>

namespace CGAL {

template < typename Base_, typename FT_, typename LA_=CGAL::LA_eigen<FT_,typename Base_::Default_ambient_dimension> >
struct Cartesian_change_FT_base : public
	Base_
{
    CGAL_CONSTEXPR Cartesian_change_FT_base(){}
    CGAL_CONSTEXPR Cartesian_change_FT_base(int d):Base_(d){}

    typedef Cartesian_change_FT_base Self;
    typedef Base_ Kernel_base;
    typedef LA_ LA;

    template <class T, class D=void> struct Type : Inherit_type<Base_, T> {};
    template <class D> struct Type <FT_tag, D> { typedef FT_ type; };
    template <class D> struct Type <RT_tag, D> { typedef FT_ type; };

    typedef NT_converter<typename Get_type<Kernel_base, FT_tag>::type,FT_> FT_converter;
    typedef transforming_iterator<FT_converter,typename Kernel_base::Point_cartesian_const_iterator> Point_cartesian_const_iterator;
    typedef transforming_iterator<FT_converter,typename Kernel_base::Vector_cartesian_const_iterator> Vector_cartesian_const_iterator;
    //FIXME: use Iterator_list!
    /*
    template<class T,bool=CGAL_BOOSTD is_same<typename iterator_tag_traits<T>::value_tag,FT_tag>::value>
    struct Iterator : Get_type<Kernel_base,T> {};
    template<class T> struct Iterator<T,true> {
      typedef transforming_iterator<FT_converter,typename Get_type<Kernel_base,T>::type> type;
    };
    */

    template<class Tag_,class Type_>
    struct Construct_cartesian_const_iterator_ {
	    typedef typename Get_functor<Kernel_base, Tag_>::type Functor_base;
	    Construct_cartesian_const_iterator_(){}
	    Construct_cartesian_const_iterator_(Self const&r):f(r){}
	    Functor_base f;
	    typedef Type_ result_type;
	    template<class T>
	    result_type operator()(T const& v, Begin_tag)const{
		    return make_transforming_iterator(f(v,Begin_tag()),FT_converter());
	    }
	    template<class T>
	    result_type operator()(T const& v, End_tag)const{
		    return make_transforming_iterator(f(v,End_tag()),FT_converter());
	    }
    };
    typedef Construct_cartesian_const_iterator_<Construct_ttag<Point_cartesian_const_iterator_tag>,Point_cartesian_const_iterator> Construct_point_cartesian_const_iterator;
    typedef Construct_cartesian_const_iterator_<Construct_ttag<Vector_cartesian_const_iterator_tag>,Vector_cartesian_const_iterator> Construct_vector_cartesian_const_iterator;

    template<class Tag_>
    struct Compute_cartesian_coordinate {
	    typedef typename Get_functor<Kernel_base, Tag_>::type Functor_base;
	    Compute_cartesian_coordinate(){}
	    Compute_cartesian_coordinate(Self const&r):f(r){}
	    Functor_base f;
	    typedef FT_ result_type;
	    template<class Obj_>
	    result_type operator()(Obj_ const& v,int i)const{
		    return FT_converter()(f(v,i));
	    }
    };

    template<class T,class U=void,class=typename Get_functor_category<Cartesian_change_FT_base,T>::type> struct Functor :
	    Inherit_functor<Kernel_base,T,U> { };
    template<class T,class U> struct Functor<T,U,Compute_tag> { };
    template<class T,class U> struct Functor<T,U,Predicate_tag> { };
    template<class D> struct Functor<Compute_point_cartesian_coordinate_tag,D,Compute_tag> {
	    typedef Compute_cartesian_coordinate<Compute_point_cartesian_coordinate_tag> type;
    };
    template<class D> struct Functor<Compute_vector_cartesian_coordinate_tag,D,Compute_tag> {
	    typedef Compute_cartesian_coordinate<Compute_vector_cartesian_coordinate_tag> type;
    };
    template<class D> struct Functor<Construct_ttag<Point_cartesian_const_iterator_tag>,D,Construct_iterator_tag> {
	    typedef Construct_point_cartesian_const_iterator type;
    };
    template<class D> struct Functor<Construct_ttag<Vector_cartesian_const_iterator_tag>,D,Construct_iterator_tag> {
	    typedef Construct_vector_cartesian_const_iterator type;
    };
};

template < typename Base_, typename FT_>
struct Cartesian_change_FT : public
	Cartesian_change_FT_base<Base_,FT_>
{
    CGAL_CONSTEXPR Cartesian_change_FT(){}
    CGAL_CONSTEXPR Cartesian_change_FT(int d):Cartesian_change_FT_base<Base_,FT_>(d){}
};

} //namespace CGAL

#endif // CGAL_KERNEL_D_CARTESIAN_CHANGE_FT_H
