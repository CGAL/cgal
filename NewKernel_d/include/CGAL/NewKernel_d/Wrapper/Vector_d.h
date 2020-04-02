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
  typedef typename Get_type<R_, RT_tag>::type                RT_;
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename R_::Kernel_base           Kbase;
  typedef typename Get_type<R_, Point_tag>::type        Point_;
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
  typedef typename Get_type<Kbase, Vector_tag>::type        Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

#if defined(BOOST_MSVC) && (BOOST_MSVC == 1900)
#  pragma warning(push)
#  pragma warning(disable: 4309)
#endif
  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Vector_d> >::value>::type> explicit Vector_d(U&&...u)
          : Rep(CVBase()(std::forward<U>(u)...)){}

#if defined(BOOST_MSVC) && (BOOST_MSVC == 1900)
#  pragma warning(pop)
#endif

//  // called from Construct_vector_d
//  template<class...U> explicit Vector_d(Eval_functor&&,U&&...u)
//          : Rep(Eval_functor(), std::forward<U>(u)...){}
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

  friend std::ostream& operator <<(std::ostream& os, const Vector_d& v)
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

  friend std::istream & operator>>(std::istream &is, Vector_d & v)
  {
    int dim;
    if( is_ascii(is) )
      is >> dim;
    else
    {
      read(is, dim);
    }
    if(!is) return is;

    std::vector<FT_> coords(dim);
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
      v = Vector_d(coords.begin(), coords.end());
    return is;
  }


  friend Vector_d operator+(const Vector_d& v,const Vector_d& w)
  {
    return typename Get_functor<R, Sum_of_vectors_tag>::type()(v,w);
  }

  friend Vector_d operator-(const Vector_d& v,const Vector_d& w)
  {
    return typename Get_functor<R, Difference_of_vectors_tag>::type()(v,w);
  }
};
#if 0
template <class R_> Vector_d<R_>::Vector_d(Vector_d &)=default;
#endif


} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_VECTOR_D_H
