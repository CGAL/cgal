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

#ifndef CGAL_WRAPPER_REF_COUNT_OBJ_H
#define CGAL_WRAPPER_REF_COUNT_OBJ_H

#include <CGAL/Origin.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#ifndef CGAL_CXX11
#include <boost/preprocessor/repetition.hpp>
#endif
#include <boost/utility/result_of.hpp>

// no need for a fancy interface here, people can use the Point_d wrapper on
// top.

namespace CGAL {

template <class R_, class Tag_>
class Ref_count_obj
{
  typedef typename R_::Kernel_base	Kbase;
  typedef typename Get_functor<Kbase, Construct_ttag<Tag_> >::type CBase;

  typedef Ref_count_obj			Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename Get_type<R_, Tag_>::type>::value));

public:
  typedef R_ R;

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  //typedef Dimension_tag<0>  Feature_dimension;

  typedef typename Get_type<Kbase, Tag_>::type	Rep;
  typedef Handle_for<Rep> Data;

private:
  Data data;
public:

  const Rep& rep() const
  {
    return CGAL::get_pointee_or_identity(data);
  }

#ifdef CGAL_CXX11
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Ref_count_obj> >::value>::type> explicit Ref_count_obj(U&&...u)
	  : data(Eval_functor(),CBase(),std::forward<U>(u)...){}

  template<class F,class...U> explicit Ref_count_obj(Eval_functor&&,F&&f,U&&...u)
	  : data(Eval_functor(),std::forward<F>(f),std::forward<U>(u)...){}

  // try not to use these
  Ref_count_obj(Rep const& v) : data(v) {}
  Ref_count_obj(Rep& v) : data(static_cast<Rep const&>(v)) {}
  Ref_count_obj(Rep&& v) : data(std::move(v)) {}

  // Do we really need this for point?
//  // this one should be implicit
//  Ref_count_obj(Origin const& v)
//    : data(Eval_functor(),CBase(),v) {}
//  Ref_count_obj(Origin& v)
//    : data(Eval_functor(),CBase(),v) {}
//  Ref_count_obj(Origin&& v)
//    : data(Eval_functor(),CBase(),std::move(v)) {}

#else

  Ref_count_obj() : data(Eval_functor(),CBase()) {}

  Ref_count_obj(Rep const& v) : data(v) {} // try not to use it

#define CGAL_CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Ref_count_obj(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : data(Eval_functor(),CBase(),BOOST_PP_ENUM_PARAMS(N,t)) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Ref_count_obj(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : data(Eval_functor(),f,BOOST_PP_ENUM_PARAMS(N,t)) {}

  BOOST_PP_REPEAT_FROM_TO(1,11,CGAL_CODE,_)
#undef CGAL_CODE
  template<class F>
  Ref_count_obj(Eval_functor,F const& f)
  : data(Eval_functor(),f) {}

//  // this one should be implicit
//  Ref_count_obj(Origin const& o)
//    : data(Eval_functor(),CBase(),o) {}

#endif

};

} //namespace CGAL

#endif
