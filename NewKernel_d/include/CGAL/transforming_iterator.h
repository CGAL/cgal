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
//
// Author(s)     : Marc Glisse

#ifndef CGAL_TRANSFORMING_ITERATOR_H
#define CGAL_TRANSFORMING_ITERATOR_H
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/utility/result_of.hpp>
#include <boost/type_traits/is_empty.hpp>
#include <CGAL/Default.h>
#include <utility>

// Inspired by the boost version, but more compact and
// without any iterator_category games.

namespace CGAL {
namespace internal {

// non-empty case
template<class T,bool=boost::is_empty<T>::value> struct Functor_as_base {
	Functor_as_base(){}
	Functor_as_base(T const& t):f(t){}
	//template<class T2> Functor_as_base(Functor_as_base<T2> const&g):f(g.functor()){}
	T const& functor()const{return f;}
	T &      functor()     {return f;}
	private:
		T f;
};

// empty case
template<class T> struct Functor_as_base<T,true> : public T {
	Functor_as_base(){}
	Functor_as_base(T const& t):T(t){}
	//template<class T2> Functor_as_base(Functor_as_base<T2> const&g):T(g.functor()){}
	T const& functor()const{return *this;}
	T &      functor()     {return *this;}
};

template <typename Derived, typename F, typename Iter, typename Ref, typename Val>
class transforming_iterator_helper
{
	typedef typename Default::Get<Ref,
#ifdef CGAL_CXX11
		decltype(std::declval<F>()(std::declval<typename std::iterator_traits<Iter>::reference>()))
#else
		typename boost::result_of<F(typename std::iterator_traits<Iter>::value_type)>::type
	// should be reference instead of value_type
#endif
			>::type reference;

	typedef typename Default::Get<Val,typename boost::remove_cv<typename boost::remove_reference<reference>::type>::type>::type value_type;

	public:
	typedef boost::iterator_adaptor<
		Derived,
		Iter,
		value_type,
		typename std::iterator_traits<Iter>::iterator_category,
		reference
			> type;
};
}

template <typename F, typename Iter, typename Ref=Default, typename Val=Default>
class transforming_iterator
: public internal::transforming_iterator_helper<transforming_iterator<F,Iter,Ref,Val>,F,Iter,Ref,Val>::type,
private internal::Functor_as_base<F>
{
	friend class boost::iterator_core_access;
	typedef typename internal::transforming_iterator_helper<transforming_iterator,F,Iter,Ref,Val>::type Base;
	typedef internal::Functor_as_base<F> Functor_base;
	typename Base::reference dereference()const{
		return functor()(*this->base_reference());
	}
	public:
	using Functor_base::functor;
	transforming_iterator(){}
	explicit transforming_iterator(Iter i,F const& f=F())
		:Base(i),Functor_base(f){}
	template<class F2,class I2,class R2,class V2>
	transforming_iterator(
		transforming_iterator<F2,I2,R2,V2> const&i,
		typename boost::enable_if_convertible<I2, Iter>::type* = 0,
		typename boost::enable_if_convertible<F2, F>::type* = 0)
		: Base(i.base()),Functor_base(i.functor()) {}

};

template <typename F, typename Iter> inline
transforming_iterator<F,Iter> make_transforming_iterator(Iter i, F const&f=F()) {
	return transforming_iterator<F,Iter>(i,f);
}

}

#endif // CGAL_TRANSFORMING_ITERATOR_H
