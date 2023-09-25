// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H 1

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/certified_quotient_predicates.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Trisegment_2.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Handle.h>
#include <CGAL/Quotient.h>
#include <CGAL/tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/Unfiltered_predicate_adaptor.h>
#include <CGAL/Lazy_exact_nt.h>

#include <boost/mpl/has_xxx.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <optional>

namespace CGAL {

namespace CGAL_SS_i {

template<class T>
T const& validate ( std::optional<T> const& o )
{
  if ( !o )
    throw std::overflow_error("Arithmetic overflow");
  return *o ;
}

template<class NT>
NT const& validate( NT const& n )
{
  if ( !CGAL_NTS is_finite(n) )
    throw std::overflow_error("Arithmetic overflow");
  return n ;
}

template<class T>
std::optional<T> cgal_make_optional( T const& v )
{
  return std::optional<T>(v) ;
}

template<class T>
std::optional<T> cgal_make_optional( bool cond, T const& v )
{
  return cond ? std::optional<T>(v) : std::optional<T>() ;
}

template<class K>
struct Is_filtering_kernel
{
  typedef Tag_false type ;
} ;

template<>
struct Is_filtering_kernel< Exact_predicates_inexact_constructions_kernel >
{
  typedef Tag_true type ;
} ;

//
// This is the same as Filtered_construction but uses optional<result> instead of exceptions.
//
template <class AC
         ,class EC
         ,class FC
         ,class C2E
         ,class C2F
         ,class E2C
         ,class F2C
         ,bool Protection = true
>
class Exceptionless_filtered_construction
{
private:
  EC Exact_construction;
  FC Filter_construction;
  C2E To_Exact;
  C2F To_Filtered;
  E2C From_Exact;
  F2C From_Filtered;

  typedef typename AC::result_type  AC_result_type;
  typedef typename FC::result_type  FC_result_type;
  typedef typename EC::result_type  EC_result_type;
  typedef typename C2F::Target_kernel FK;

  bool has_enough_precision(const typename FK::Point_2& point, double precision) const
  {
    return has_smaller_relative_precision(point.x(), precision) &&
           has_smaller_relative_precision(point.y(), precision);
  }

  bool has_enough_precision(const std::tuple<typename FK::FT, typename FK::Point_2>& time_and_point, double precision) const
  {
    return has_smaller_relative_precision(std::get<0>(time_and_point), precision) &&
           has_enough_precision(std::get<1>(time_and_point), precision);
  }

  bool has_enough_precision(const CGAL::Trisegment_2<FK, CGAL_SS_i::Segment_2_with_ID<FK> >& trisegment, double precision) const
  {
    return has_enough_precision(trisegment.e0().source(), precision) &&
           has_enough_precision(trisegment.e0().target(), precision) &&
           has_smaller_relative_precision(trisegment.w0(), precision) &&
           has_enough_precision(trisegment.e1().source(), precision) &&
           has_enough_precision(trisegment.e1().target(), precision) &&
           has_smaller_relative_precision(trisegment.w1(), precision) &&
           has_enough_precision(trisegment.e2().source(), precision) &&
           has_enough_precision(trisegment.e2().target(), precision) &&
           has_smaller_relative_precision(trisegment.w2(), precision);
  }

public:

  Exceptionless_filtered_construction() {}

  Exceptionless_filtered_construction(const EC& Exact_construction, const FC& Filter_construction)
    : Exact_construction(Exact_construction)
    , Filter_construction(Filter_construction)
  {}

  typedef AC_result_type           result_type;

  template <class ... A>
  result_type
  operator()(A&& ... a) const
  {
    {
      Protect_FPU_rounding<Protection> P;
      try
      {
        FC_result_type fr = Filter_construction(To_Filtered(std::forward<A>(a))...);

        const double precision =
          Lazy_exact_nt<double>::get_relative_precision_of_to_double();

        if ( fr && has_enough_precision(*fr, precision) )
          return From_Filtered(fr);
      }
      catch (Uncertain_conversion_exception&) {}
    }

    Protect_FPU_rounding<!Protection> P(CGAL_FE_TONEAREST);
    CGAL_expensive_assertion(FPU_get_cw() == CGAL_FE_TONEAREST);
    EC_result_type er = Exact_construction(To_Exact(std::forward<A>(a))...) ;
    return From_Exact(er);
  }
};


//
// This number type is provided because unlike Quotient<> it allows you to create it
// with a zero denominator. Of course you can't evaluate it in that case, but is convenient because it allows client code
// to handle the "error" itself, which in this context is useful.
//
template<class NT>
class Rational
{
  public:

    Rational( NT aN, NT aD ) : mN(aN), mD(aD) {}

    NT n() const { return mN ; }
    NT d() const { return mD ; }

    CGAL::Quotient<NT> to_quotient() const { return CGAL::Quotient<NT>(mN,mD) ; }

    NT to_nt() const { return mN / mD ; }

    friend std::ostream& operator << ( std::ostream& os, Rational<NT> const& rat )
    {
      if ( ! CGAL_NTS is_zero(rat.d()) )
           return os << n2str(rat.n()/rat.d());
      else return os << "INF_RATIONAL" ;
    }

  private:

    NT mN, mD ;
} ;


// Debug struct
template <class Info>
struct FPU_checker;

template <class K>
struct FPU_checker<std::optional< Line_2<K> > >
{
  static bool is_valid()
  {
    return !std::is_same<typename K::FT, CGAL::Interval_nt<false> >::value ||
           FPU_get_cw() == CGAL_FE_UPWARD;
  }
};

template <class FT>
struct FPU_checker<std::optional< CGAL_SS_i::Rational< FT > > >
{
  static bool is_valid()
  {
    return !std::is_same<FT, CGAL::Interval_nt<false> >::value ||
           FPU_get_cw() == CGAL_FE_UPWARD;
  }
};

template <class K>
struct FPU_checker<std::optional< Point_2<K> > >
{
  static bool is_valid()
  {
    return !std::is_same<typename K::FT, CGAL::Interval_nt<false> >::value ||
           FPU_get_cw() == CGAL_FE_UPWARD;
  }
};


template<class K>
struct Functor_base_2
{
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  typedef typename K::Vector_2 Vector_2 ;
  typedef CGAL_SS_i::Segment_2_with_ID<K> Segment_2_with_ID ;

  typedef CGAL::Trisegment_2<K, Segment_2_with_ID> Trisegment_2 ;

  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;
};

template<class Converter>
struct SS_converter : Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Source_kernel::FT Source_FT ;
  typedef typename Target_kernel::FT Target_FT ;

  typedef typename Source_kernel::Point_2 Source_point_2 ;
  typedef typename Target_kernel::Point_2 Target_point_2 ;

  typedef typename Source_kernel::Vector_2 Source_vector_2 ;
  typedef typename Target_kernel::Vector_2 Target_vector_2 ;

  typedef typename Source_kernel::Segment_2 Source_segment_2 ;
  typedef typename Target_kernel::Segment_2 Target_segment_2 ;

  typedef Segment_2_with_ID<Source_kernel> Source_segment_2_with_ID ;
  typedef Segment_2_with_ID<Target_kernel> Target_segment_2_with_ID ;

  typedef Trisegment_2<Source_kernel, Source_segment_2_with_ID> Source_trisegment_2 ;
  typedef Trisegment_2<Target_kernel, Target_segment_2_with_ID> Target_trisegment_2 ;

  typedef std::tuple<Source_FT,Source_point_2> Source_time_and_point_2 ;
  typedef std::tuple<Target_FT,Target_point_2> Target_time_and_point_2 ;

  typedef std::optional<Source_FT> Source_opt_FT ;
  typedef std::optional<Target_FT> Target_opt_FT ;

  typedef std::optional<Source_point_2> Source_opt_point_2 ;
  typedef std::optional<Target_point_2> Target_opt_point_2 ;

  typedef std::optional<Source_time_and_point_2> Source_opt_time_and_point_2 ;
  typedef std::optional<Target_time_and_point_2> Target_opt_time_and_point_2 ;

  typedef std::optional<Source_segment_2> Source_opt_segment_2 ;
  typedef std::optional<Target_segment_2> Target_opt_segment_2 ;

  typedef typename Source_trisegment_2::Self_ptr Source_trisegment_2_ptr ;
  typedef typename Target_trisegment_2::Self_ptr Target_trisegment_2_ptr ;


  Target_FT cvt_n(Source_FT const& n) const  { return this->Converter::operator()(n); }

  Target_opt_FT cvt_n(Source_opt_FT const& n) const
  {
    Target_opt_FT r ;
    if ( n )
      r = cvt_n(*n);
    return r ;
  }

  Target_point_2   cvt_p(Source_point_2 const& p) const  { return this->Converter::operator()(p); }

  Target_vector_2 cvt_v( Source_vector_2 const& v) const {
    return Target_vector_2(cvt_p(Source_point_2(CGAL::ORIGIN)), cvt_p(Source_point_2(CGAL::ORIGIN) + v) ) ;
  }

  Target_segment_2 cvt_s( Source_segment_2 const& e) const {
    return Target_segment_2(cvt_p(e.source()), cvt_p(e.target())) ;
  }

  Target_segment_2_with_ID cvt_s( Source_segment_2_with_ID const& e) const {
    return Target_segment_2_with_ID(cvt_p(e.source()), cvt_p(e.target()), e.mID) ;
  }

  Target_time_and_point_2 cvt_t_p( Source_time_and_point_2 const& v ) const
  {
    Source_FT      t ;
    Source_point_2 p ;
    std::tie(t,p) = v ;
    return Target_time_and_point_2(cvt_n(t),cvt_p(p));
  }

  Target_trisegment_2_ptr cvt_single_trisegment( Source_trisegment_2_ptr const& tri ) const
  {
    CGAL_precondition( tri!= Source_trisegment_2_ptr() ) ;

    return Target_trisegment_2_ptr ( new Target_trisegment_2(cvt_s(tri->e0())
                                                            ,cvt_n(tri->w0())
                                                            ,cvt_s(tri->e1())
                                                            ,cvt_n(tri->w1())
                                                            ,cvt_s(tri->e2())
                                                            ,cvt_n(tri->w2())
                                                            ,tri->collinearity()
                                                            ,tri->id()
                                                            )
                                   ) ;
  }

  Target_trisegment_2_ptr cvt_trisegment( Source_trisegment_2_ptr const& tri ) const
  {
    Target_trisegment_2_ptr res ;

    if ( tri )
    {
      res = cvt_single_trisegment(tri) ;

      if ( tri->child_l() )
        res->set_child_l( cvt_trisegment(tri->child_l()) ) ;

      if ( tri->child_r() )
        res->set_child_r( cvt_trisegment(tri->child_r() ) ) ;

      if ( tri->child_t() )
        res->set_child_t( cvt_trisegment(tri->child_t() ) ) ;
    }

    return res ;
  }

  bool operator()( bool v ) const { return v ; }

  Trisegment_collinearity  operator()(Trisegment_collinearity c) const { return c ; }

  Oriented_side operator()(Oriented_side s) const { return s ; }

  Target_FT        operator()(Source_FT const& n) const { return cvt_n(n) ; }

  Target_opt_FT    operator()(Source_opt_FT const& n) const { return cvt_n(n) ; }

  Target_point_2   operator()( Source_point_2 const& p) const { return cvt_p(p) ; }

  Target_vector_2 operator()( Source_vector_2 const& v) const { return cvt_v(v); }

  Target_segment_2 operator()( Source_segment_2 const& s) const { return cvt_s(s); }

  Target_segment_2_with_ID operator()( Source_segment_2_with_ID const& s) const { return cvt_s(s); }

  Target_trisegment_2_ptr operator()( Source_trisegment_2_ptr const& tri ) const
  {
    return cvt_trisegment(tri);
  }

  Target_time_and_point_2 operator() ( Source_time_and_point_2 const& v ) const
  {
    return cvt_t_p(v);
  }

  Target_opt_point_2 operator()( Source_opt_point_2 const& p) const
  {
    if ( p )
         return Target_opt_point_2(cvt_p(*p));
    else return Target_opt_point_2();
  }

  Target_opt_segment_2 operator()( Source_opt_segment_2 const& s) const
  {
    if ( s )
         return Target_opt_segment_2(cvt_s(*s));
    else return Target_opt_segment_2();
  }

  Target_opt_time_and_point_2 operator()( Source_opt_time_and_point_2 const& v) const
  {
    if ( v )
         return Target_opt_time_and_point_2(cvt_t_p(*v));
    else return Target_opt_time_and_point_2();
  }

};

BOOST_MPL_HAS_XXX_TRAIT_DEF(Filters_split_events_tag)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Protector)
BOOST_MPL_HAS_XXX_TRAIT_DEF(Segment_2_with_ID)

template <class GT, bool has_Protector = has_Protector<GT>::value>
struct Get_protector{ struct type{}; };

template <class GT>
struct Get_protector<GT, true>
{
  typedef typename GT::Protector type;
};

} // namespace CGAL_SS_i
} // namespace CGAL

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H
