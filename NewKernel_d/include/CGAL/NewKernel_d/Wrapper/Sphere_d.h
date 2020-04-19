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

#ifndef CGAL_WRAPPER_SPHERE_D_H
#define CGAL_WRAPPER_SPHERE_D_H

#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Wrap {

template <class R_>
class Sphere_d : public Get_type<typename R_::Kernel_base, Sphere_tag>::type
{
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename R_::Kernel_base                Kbase;
  typedef typename Get_type<R_, Point_tag>::type        Point_;
  typedef typename Get_functor<Kbase, Construct_ttag<Sphere_tag> >::type        CSBase;
  typedef typename Get_functor<Kbase, Center_of_sphere_tag>::type                COSBase;
  typedef typename Get_functor<Kbase, Squared_radius_tag>::type                        SRBase;

  typedef Sphere_d                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename Get_type<R_, Sphere_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef typename Increment_dimension<Ambient_dimension,-1>::type Feature_dimension;

  typedef typename Get_type<Kbase, Sphere_tag>::type        Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Sphere_d> >::value>::type> explicit Sphere_d(U&&...u)
          : Rep(CSBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//          : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Sphere_d(Eval_functor&&,F&&f,U&&...u)
          : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Sphere_d(Rep const& v) : Rep(v) {}
  Sphere_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Sphere_d(Rep&& v) : Rep(std::move(v)) {}


    //TODO: if COSBase returns a reference to a base point, cast it to a
    //reference to a wrapper point. Ugly but should be safe.
    Point_ center()const{
      return Point_(Eval_functor(),COSBase(),rep());
    }
  FT_ squared_radius()const{
    return SRBase()(rep());
  }

};

} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_SPHERE_D_H
