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

#ifndef CGAL_STORE_KERNEL_H
#define CGAL_STORE_KERNEL_H

#include <CGAL/assertions.h>
#include <boost/type_traits/is_empty.hpp>

namespace CGAL {
namespace internal {
BOOST_MPL_HAS_XXX_TRAIT_DEF(Do_not_store_kernel)
template<class T,bool=boost::is_empty<T>::value,bool=has_Do_not_store_kernel<T>::value> struct Do_not_store_kernel {
	enum { value=false };
	typedef Tag_false type;
};
template<class T> struct Do_not_store_kernel<T,true,false> {
	enum { value=true };
	typedef Tag_true type;
};
template<class T,bool b> struct Do_not_store_kernel<T,b,true> {
	typedef typename T::Do_not_store_kernel type;
	enum { value=type::value };
};
}

template<class R_,bool=internal::Do_not_store_kernel<R_>::value>
struct Store_kernel {
	Store_kernel(){}
	Store_kernel(R_ const&){}
	enum { kernel_is_stored = false };
	R_ kernel()const{return R_();}
	typedef R_ reference_type;
	void set_kernel(R_ const&){}
};
template<class R_>
struct Store_kernel<R_,false> {
	Store_kernel():rp(0){
		CGAL_warning_msg(true,"I should know my kernel");
	}
	Store_kernel(R_ const& r):rp(&r){}
	enum { kernel_is_stored = true };
	R_ const& kernel()const{
		CGAL_warning_msg(rp!=0,"I should know my kernel");
		return *rp;
	}
	typedef R_ const& reference_type;
	void set_kernel(R_ const&r){rp=&r;}
	private:
	R_ const* rp;
};

//For a second kernel. TODO: find something more elegant
template<class R_,bool=internal::Do_not_store_kernel<R_>::value>
struct Store_kernel2 {
	Store_kernel2(){}
	Store_kernel2(R_ const&){}
	enum { kernel2_is_stored = false };
	R_ kernel2()const{return R_();}
	typedef R_ reference2_type;
	void set_kernel2(R_ const&){}
};
template<class R_>
struct Store_kernel2<R_,false> {
	Store_kernel2(){
		//CGAL_warning_msg(true,"I should know my kernel");
	}
	Store_kernel2(R_ const& r):rp(&r){}
	enum { kernel2_is_stored = true };
	R_ const& kernel2()const{
		CGAL_warning_msg(rp==0,"I should know my kernel");
		return *rp;
	}
	typedef R_ const& reference2_type;
	void set_kernel2(R_ const&r){rp=&r;}
	private:
	R_ const* rp;
};
}
#define CGAL_BASE_INIT(X,Y) \
	X():Y(){} \
	X(R_ const&r):Y(r){}
#define CGAL_FUNCTOR_INIT_STORE(X) CGAL_BASE_INIT(X,Store_kernel<R_>)
#define CGAL_FUNCTOR_INIT_IGNORE(X) \
	X(){} \
	X(R_ const&){}

#endif // CGAL_STORE_KERNEL_H
