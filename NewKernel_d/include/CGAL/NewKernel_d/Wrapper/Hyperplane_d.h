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

#ifndef CGAL_WRAPPER_HYPERPLANE_D_H
#define CGAL_WRAPPER_HYPERPLANE_D_H

#include <CGAL/representation_tags.h>
#include <CGAL/assertions.h>
#include <boost/type_traits.hpp>
#include <CGAL/Kernel/Return_base_tag.h>
#include <CGAL/Dimension.h>

namespace CGAL {
namespace Wrap {

template <class R_>
class Hyperplane_d : public Get_type<typename R_::Kernel_base, Hyperplane_tag>::type
{
  typedef typename Get_type<R_, FT_tag>::type                FT_;
  typedef typename R_::Kernel_base                Kbase;
  typedef typename Get_type<R_, Vector_tag>::type        Vector_;
  typedef typename Get_functor<Kbase, Construct_ttag<Hyperplane_tag> >::type        CHBase;
  typedef typename Get_functor<Kbase, Orthogonal_vector_tag>::type                OVBase;
  typedef typename Get_functor<Kbase, Hyperplane_translation_tag>::type                        HTBase;

  typedef Hyperplane_d                            Self;
  CGAL_static_assertion((boost::is_same<Self, typename Get_type<R_, Hyperplane_tag>::type>::value));

public:

  typedef Tag_true Is_wrapper;
  typedef typename R_::Default_ambient_dimension Ambient_dimension;
  typedef typename Increment_dimension<Ambient_dimension,-1>::type Feature_dimension;

  typedef typename Get_type<Kbase, Hyperplane_tag>::type        Rep;

  const Rep& rep() const
  {
    return *this;
  }

  Rep& rep()
  {
    return *this;
  }

  typedef          R_                       R;

  template<class...U,class=typename std::enable_if<!std::is_same<std::tuple<typename std::decay<U>::type...>,std::tuple<Hyperplane_d> >::value>::type> explicit Hyperplane_d(U&&...u)
          : Rep(CHBase()(std::forward<U>(u)...)){}

//  // called from Construct_point_d
//  template<class...U> explicit Point_d(Eval_functor&&,U&&...u)
//          : Rep(Eval_functor(), std::forward<U>(u)...){}
  template<class F,class...U> explicit Hyperplane_d(Eval_functor&&,F&&f,U&&...u)
          : Rep(std::forward<F>(f)(std::forward<U>(u)...)){}

#if 0
  // the new standard may make this necessary
  Point_d(Point_d const&)=default;
  Point_d(Point_d &);//=default;
  Point_d(Point_d &&)=default;
#endif

  // try not to use these
  Hyperplane_d(Rep const& v) : Rep(v) {}
  Hyperplane_d(Rep& v) : Rep(static_cast<Rep const&>(v)) {}
  Hyperplane_d(Rep&& v) : Rep(std::move(v)) {}


  //TODO: if OVBase returns a reference to a base vector, cast it to a
  //reference to a wrapper vector. Ugly but should be safe.
  Vector_ orthogonal_vector()const{
    return Vector_(Eval_functor(),OVBase(),rep());
  }
  FT_ translation()const{
    return HTBase()(rep());
  }


};

} //namespace Wrap
} //namespace CGAL

#endif // CGAL_WRAPPER_SPHERE_D_H
