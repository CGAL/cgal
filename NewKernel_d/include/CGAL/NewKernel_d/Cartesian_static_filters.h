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

#ifndef CGAL_KD_CARTESIAN_STATIC_FILTERS_H
#define CGAL_KD_CARTESIAN_STATIC_FILTERS_H
#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/Dimension.h>
#include <CGAL/internal/Static_filters/tools.h> // bug, should be included by the next one
#include <CGAL/internal/Static_filters/Orientation_2.h>
#include <CGAL/internal/Static_filters/Side_of_oriented_circle_2.h>
#include <boost/mpl/if.hpp>

namespace CGAL {
namespace SFA { // static filter adapter
// Note that this would be quite a bit simpler without stateful kernels
template <class Base_,class R_> struct Adapter_2 {
        typedef typename Get_type<R_, Point_tag>::type        Point;
        typedef typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type CC;
        typedef typename Get_functor<Base_, Orientation_of_points_tag>::type Orientation_base;
        typedef typename Get_functor<Base_, Side_of_oriented_sphere_tag>::type Side_of_oriented_circle_base;
        struct Point_2 {
                R_ const&r; CC const&c; Point const& p;
                Point_2(R_ const&r_, CC const&c_, Point const&p_):r(r_),c(c_),p(p_){}
                decltype(auto) x()const{return c(p,0);}
                decltype(auto) y()const{return c(p,1);}
        };
        struct Vector_2 {};
        struct Circle_2 {};
        struct Orientation_2 {
                typedef typename Get_type<R_, Orientation_tag>::type result_type;
                auto operator()(Point_2 const&A, Point_2 const&B, Point_2 const&C)const{
                        Point const* t[3]={&A.p,&B.p,&C.p};
                        return Orientation_base(A.r)(make_transforming_iterator<Dereference_functor>(t+0),make_transforming_iterator<Dereference_functor>(t+3));
                }
        };
        struct Side_of_oriented_circle_2 {
                typedef typename Get_type<R_, Oriented_side_tag>::type result_type;
                auto operator()(Point_2 const&A, Point_2 const&B, Point_2 const&C, Point_2 const&D)const{
                        Point const* t[3]={&A.p,&B.p,&C.p};
                        return Side_of_oriented_circle_base(A.r)(make_transforming_iterator<Dereference_functor>(t+0),make_transforming_iterator<Dereference_functor>(t+3),D.p);
                }
        };
};
template <class Base_,class R_> struct Orientation_of_points_2 : private Store_kernel<R_> {
        CGAL_FUNCTOR_INIT_STORE(Orientation_of_points_2)
        typedef typename Get_type<R_, Point_tag>::type        Point;
        typedef typename Get_type<R_, Orientation_tag>::type result_type;
        typedef typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type CC;
        typedef Adapter_2<Base_, R_> Adapter;
        template<class Iter> result_type operator()(Iter f, Iter CGAL_assertion_code(e))const{
                CC c(this->kernel());
                Point const& A=*f;
                Point const& B=*++f;
                Point const& C=*++f;
                CGAL_assertion(++f==e);
                typedef typename Adapter::Point_2 P;
                return typename internal::Static_filters_predicates::Orientation_2<Adapter>()(P(this->kernel(),c,A),P(this->kernel(),c,B),P(this->kernel(),c,C));
        }
};
template <class Base_,class R_> struct Side_of_oriented_sphere_2 : private Store_kernel<R_> {
        CGAL_FUNCTOR_INIT_STORE(Side_of_oriented_sphere_2)
        typedef typename Get_type<R_, Point_tag>::type        Point;
        typedef typename Get_type<R_, Oriented_side_tag>::type result_type;
        typedef typename Get_functor<R_, Compute_point_cartesian_coordinate_tag>::type CC;
        typedef Adapter_2<Base_, R_> Adapter;
        template<class Iter> result_type operator()(Iter f, Iter CGAL_assertion_code(e), Point const& D)const{
                CC c(this->kernel());
                Point const& A=*f;
                Point const& B=*++f;
                Point const& C=*++f;
                CGAL_assertion(++f==e);
                typedef typename Adapter::Point_2 P;
                return typename internal::Static_filters_predicates::Side_of_oriented_circle_2<Adapter>()(P(this->kernel(),c,A),P(this->kernel(),c,B),P(this->kernel(),c,C),P(this->kernel(),c,D));
        }
};
}

template <class Dim_ /* should be implicit */, class R_, class Derived_=Default>
struct Cartesian_static_filters : public R_ {
  constexpr Cartesian_static_filters(){}
  constexpr Cartesian_static_filters(int d):R_(d){}
};

template <class R_, class Derived_>
struct Cartesian_static_filters<Dimension_tag<2>, R_, Derived_> : public R_ {
  constexpr Cartesian_static_filters(){}
  constexpr Cartesian_static_filters(int d):R_(d){}
        typedef Cartesian_static_filters<Dimension_tag<2>, R_, Derived_> Self;
        typedef typename Default::Get<Derived_,Self>::type Derived;
        template <class T, class=void> struct Functor : Inherit_functor<R_, T> {};
        template <class D> struct Functor <Orientation_of_points_tag,D> {
                typedef
                        //typename boost::mpl::if_ <
                        //boost::is_same<D,No_filter_tag>,
                        //typename Get_functor<R_, Orientation_of_points_tag>::type,
                        SFA::Orientation_of_points_2<R_,Derived>
                        //        >::type
                                type;
        };
        template <class D> struct Functor <Side_of_oriented_sphere_tag,D> {
                typedef SFA::Side_of_oriented_sphere_2<R_,Derived> type;
        };
};

}

#endif
