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

#ifndef CGAL_WRAPPER_WEIGHTED_POINT_D_H
#define CGAL_WRAPPER_WEIGHTED_POINT_D_H

#include <istream>
#include <ostream>
#include <CGAL/IO/io.h>
#include <CGAL/representation_tags.h>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Wrap {

template <class R_>
class Weighted_point_d : public Get_type<typename R_::Kernel_base, Weighted_point_tag>::type
{
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename R_::Kernel_base                Kbase;
  typedef typename Get_type<R_, Point_tag>::type        Point_;
  typedef typename Get_functor<Kbase, Construct_ttag<Weighted_point_tag> >::type        CWPBase;
  typedef typename Get_functor<Kbase, Point_drop_weight_tag>::type                PDWBase;
  typedef typename Get_functor<Kbase, Point_weight_tag>::type                        PWBase;

  typedef Weighted_point_d                            Self;
  BOOST_STATIC_ASSERT((boost::is_same<Self, typename Get_type<R_, Weighted_point_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef Dimension_tag<0> Feature_dimension;

  typedef typename Get_type<Kbase, Weighted_point_tag>::type        Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Weighted_point_d> >::value>::type> explicit Weighted_point_d(U&&...u)
          : Rep(CWPBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//          : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Weighted_point_d(Eval_functor&&,F&&f,U&&...u)
          : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Weighted_point_d(Rep const& v) : Rep(v) {}
  Weighted_point_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Weighted_point_d(Rep&& v) : Rep(std::move(v)) {}


  //TODO: use references?
  Point_ point()const{
    return Point_(Eval_functor(),PDWBase(),rep());
  }
  FT_ weight()const{
    return PWBase()(rep());
  }

};

template <class R_>
std::ostream& operator<<(std::ostream& os, const Weighted_point_d<R_>& p)
{
  if(is_ascii(os))
  {
    return os << p.point() << ' ' << p.weight();
  }
  else
  {
    write(os, p.point());
    write(os, p.weight());
  }
  return os;
}

template <class R_>
std::istream& operator>>(std::istream& is, Weighted_point_d<R_>& p)
{
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename Get_type<R_, Point_tag>::type        Point_;
  typedef typename Get_functor<R_, Construct_ttag<Weighted_point_tag> >::type        CWP;
  Point_ q; FT_ w;
  if(is_ascii(is))
  {
    if(is >> q >> iformat(w)) p=CWP()(q,w);
  }
  else
  {
    read(is, q);
    read(is, w);
    p=CWP()(q,w);
  }
  return is;
}

} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_WEIGHTED_POINT_D_H
