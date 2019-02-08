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

#ifndef CGAL_WRAPPER_VECTOR_D_H
#define CGAL_WRAPPER_VECTOR_D_H

#include <istream>
#include <ostream>
#include <CGAL/IO/io.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/mpl.h>
#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>
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


  decltype(auto) cartesian(int i)const{
	  return CCBase()(rep(),i);
  }

  decltype(auto) operator[](int i)const{
	  return CCBase()(rep(),i);
  }

  decltype(auto) cartesian_begin()const{
	  return CVI()(rep(),Begin_tag());
  }

  decltype(auto) cartesian_end()const{
	  return CVI()(rep(),End_tag());
  }

  Vector_d operator-() const
  {
    return typename Get_functor<R, Opposite_vector_tag>::type()(*this);
  }

  int dimension() const {
    typedef typename Get_functor<Kbase, Vector_dimension_tag>::type VDBase;
    return VDBase()(rep());
  }

  decltype(auto) squared_length()const{
	  return SLBase()(rep());
  }
};
#if 0
template <class R_> Vector_d<R_>::Vector_d(Vector_d &)=default;
#endif

template <class R_>
std::ostream& operator <<(std::ostream& os, const Vector_d<R_>& v)
{
  auto b = v.cartesian_begin();
  auto e = v.cartesian_end();
  if(is_ascii(os))
  {
    os << v.dimension();
    for(; b != e; ++b){
      os << " " << *b;
    }
  }
  else
  {
    write(os, v.dimension());
    for(; b != e; ++b){
      write(os, *b);
    }
  }
  
  return os;
}


template<typename K>
std::istream &
operator>>(std::istream &is, Vector_d<K> & v)
{
  typedef typename Get_type<K, Vector_tag>::type V;
  typedef typename Get_type<K, FT_tag>::type   FT;
  int dim;
  if( is_ascii(is) )
    is >> dim;
  else
  {
    read(is, dim);
  }
  if(!is) return is;

  std::vector<FT> coords(dim);
  if(is_ascii(is))
  {
    for(int i=0;i<dim;++i)
      is >> iformat(coords[i]);
  }
  else
  {
    for(int i=0;i<dim;++i)
      read(is, coords[i]);
  }

  if(is)
    v = V(coords.begin(), coords.end());
  return is;
}


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
