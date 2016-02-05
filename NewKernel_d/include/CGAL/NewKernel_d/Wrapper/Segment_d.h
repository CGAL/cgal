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

#ifndef CGAL_WRAPPER_SEGMENT_D_H
#define CGAL_WRAPPER_SEGMENT_D_H

#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
#ifndef CGAL_CXX11
#include <boost/preprocessor/repetition.hpp>
#endif
#include <boost/utility/result_of.hpp>

namespace CGAL {
namespace Wrap {

template <class R_>
class Segment_d : public Get_type<typename R_::Kernel_base, Segment_tag>::type
{
  typedef typename Get_type<R_, RT_tag>::type		RT_;
  typedef typename Get_type<R_, FT_tag>::type		FT_;
  typedef typename R_::Kernel_base			Kbase;
  typedef typename Get_type<R_, Point_tag>::type	Point_;
  typedef typename Get_functor<Kbase, Construct_ttag<Point_tag> >::type CPBase;
  typedef typename Get_functor<Kbase, Construct_ttag<Segment_tag> >::type CSBase;
  typedef typename Get_functor<Kbase, Segment_extremity_tag>::type CSEBase;

  typedef Segment_d                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename Get_type<R_, Segment_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<1>  Feature_dimension;

  typedef typename Get_type<Kbase, Segment_tag>::type	Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

#ifdef CGAL_CXX11
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Segment_d> >::value>::type> explicit Segment_d(U&&...u)
	  : Rep(CSBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//	  : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Segment_d(Eval_functor&&,F&&f,U&&...u)
	  : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Segment_d(Rep const& v) : Rep(v) {}
  Segment_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Segment_d(Rep&& v) : Rep(std::move(v)) {}

#else

  Segment_d() : Rep(CSBase()()) {}

  Segment_d(Rep const& v) : Rep(v) {} // try not to use it

#define CGAL_CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Segment_d(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(CSBase()( \
	BOOST_PP_ENUM_PARAMS(N,t))) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Segment_d(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(f(BOOST_PP_ENUM_PARAMS(N,t))) {}
  /*
  template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  Point_d(Eval_functor,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(Eval_functor(), BOOST_PP_ENUM_PARAMS(N,t)) {}
  */

  BOOST_PP_REPEAT_FROM_TO(1,11,CGAL_CODE,_)
#undef CGAL_CODE

#endif

	  //TODO: if CSEBase returns a reference to a base point, cast it to a
	  //reference to a wrapper point. Ugly but should be safe.
	  Point_ source()const{
		  return Point_(Eval_functor(),CSEBase(),rep(),0);
	  }
	  Point_ target()const{
		  return Point_(Eval_functor(),CSEBase(),rep(),1);
	  }

};

} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_SEGMENT_D_H
