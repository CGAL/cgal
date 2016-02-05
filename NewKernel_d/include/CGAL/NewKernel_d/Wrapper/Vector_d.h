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

#ifndef CGAL_WRAPPER_VECTOR_D_H
#define CGAL_WRAPPER_VECTOR_D_H

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
class Vector_d : public Get_type<typename R_::Kernel_base, Vector_tag>::type
{
  typedef typename Get_type<R_, RT_tag>::type		RT_;
  typedef typename Get_type<R_, FT_tag>::type		FT_;
  typedef typename R_::Kernel_base           Kbase;
  typedef typename Get_type<R_, Point_tag>::type	Point_;
  typedef typename Get_functor<Kbase, Construct_ttag<Vector_tag> >::type CVBase;
  typedef typename Get_functor<Kbase, Compute_vector_cartesian_coordinate_tag>::type CCBase;
  typedef typename Get_functor<Kbase, Construct_ttag<Vector_cartesian_const_iterator_tag> >::type CVI;
  typedef typename Get_functor<Kbase, Squared_length_tag>::type SLBase;

  typedef Vector_d                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename Get_type<R_, Vector_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<0>  Feature_dimension;

  //typedef typename R_::Vector_cartesian_const_iterator Cartesian_const_iterator;
  typedef typename Get_type<Kbase, Vector_tag>::type	Rep;

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
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Vector_d> >::value>::type> explicit Vector_d(U&&...u)
	  : Rep(CVBase()(std::forward<U>(u)...)){}

//  // called from Construct_vector_d
//  template<class...U> explicit Vector_d(Eval_functor&&,U&&...u)
//	  : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Vector_d(Eval_functor&&,F&&f,U&&...u)
	  : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Vector_d(Vector_d const&)=default;
  Vector_d(Vector_d &);//=default;
  Vector_d(Vector_d &&)=default;
#endif

  // try not to use these
  Vector_d(Rep const& v) : Rep(v) {}
  Vector_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Vector_d(Rep&& v) : Rep(std::move(v)) {}

  // this one should be implicit
  Vector_d(Null_vector const& v)
    : Rep(CVBase()(v)) {}
  Vector_d(Null_vector& v)
    : Rep(CVBase()(v)) {}
  Vector_d(Null_vector&& v)
    : Rep(CVBase()(std::move(v))) {}

#else

  Vector_d() : Rep(CVBase()()) {}

  Vector_d(Rep const& v) : Rep(v) {} // try not to use it

#define CGAL_CODE(Z,N,_) template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  explicit Vector_d(BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(CVBase()( \
	BOOST_PP_ENUM_PARAMS(N,t))) {} \
  \
  template<class F,BOOST_PP_ENUM_PARAMS(N,class T)> \
  Vector_d(Eval_functor,F const& f,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(f(BOOST_PP_ENUM_PARAMS(N,t))) {}
  /*
  template<BOOST_PP_ENUM_PARAMS(N,class T)> \
  Vector_d(Eval_functor,BOOST_PP_ENUM_BINARY_PARAMS(N,T,const&t)) \
  : Rep(Eval_functor(), BOOST_PP_ENUM_PARAMS(N,t)) {}
  */

  BOOST_PP_REPEAT_FROM_TO(1,11,CGAL_CODE,_)
#undef CGAL_CODE

  // this one should be implicit
  Vector_d(Null_vector const& v)
    : Rep(CVBase()(v)) {}

#endif

  typename boost::result_of<CCBase(Rep,int)>::type cartesian(int i)const{
	  return CCBase()(rep(),i);
  }

  typename boost::result_of<CCBase(Rep,int)>::type operator[](int i)const{
	  return CCBase()(rep(),i);
  }

  typename boost::result_of<CVI(Rep,Begin_tag)>::type cartesian_begin()const{
	  return CVI()(rep(),Begin_tag());
  }

  typename boost::result_of<CVI(Rep,End_tag)>::type cartesian_end()const{
	  return CVI()(rep(),End_tag());
  }

  Vector_d operator-() const
  {
    return typename Get_functor<R, Opposite_vector_tag>::type()(*this);
  }

  /*
  Direction_d direction() const
  {
    return R().construct_direction_d_object()(*this);
  }

  Vector_d transform(const Aff_transformation_d &t) const
  {
    return t.transform(*this);
  }

  Vector_d operator/(const RT& c) const
  {
   return R().construct_divided_vector_d_object()(*this,c);
  }

  Vector_d operator/(const typename First_if_different<FT_,RT>::Type & c) const
  {
   return R().construct_divided_vector_d_object()(*this,c);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  x() const
  {
    return R().compute_x_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_y_3, Vector_3>::type
  y() const
  {
    return R().compute_y_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_z_3, Vector_3>::type
  z() const
  {
    return R().compute_z_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hx_3, Vector_3>::type
  hx() const
  {
    return R().compute_hx_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hy_3, Vector_3>::type
  hy() const
  {
    return R().compute_hy_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hz_3, Vector_3>::type
  hz() const
  {
    return R().compute_hz_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
  hw() const
  {
    return R().compute_hw_3_object()(*this);
  }

  typename Qualified_result_of<typename R::Compute_x_3, Vector_3>::type
  cartesian(int i) const
  {
    CGAL_kernel_precondition( (i == 0) || (i == 1) || (i == 2) );
    if (i==0) return x();
    if (i==1) return y();
    return z();
  }

  typename Qualified_result_of<typename R::Compute_hw_3, Vector_3>::type
  homogeneous(int i) const
  {
    CGAL_kernel_precondition( (i >= 0) || (i <= 3) );
    if (i==0) return hx();
    if (i==1) return hy();
    if (i==2) return hz();
    return hw();
  }

  int dimension() const // bad idea?
  {
      return rep.dimension();
  }
*/
  typename boost::result_of<SLBase(Rep)>::type squared_length()const{
	  return SLBase()(rep());
  }
};
#if 0
template <class R_> Vector_d<R_>::Vector_d(Vector_d &)=default;
#endif

//TODO: IO

template <class R_>
Vector_d<R_> operator+(const Vector_d<R_>& v,const Vector_d<R_>& w)
{
	return typename Get_functor<R_, Sum_of_vectors_tag>::type()(v,w);
}

template <class R_>
Vector_d<R_> operator-(const Vector_d<R_>& v,const Vector_d<R_>& w)
{
	return typename Get_functor<R_, Difference_of_vectors_tag>::type()(v,w);
}

} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_VECTOR_D_H
