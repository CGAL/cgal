// Copyright (c) 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H 1

#include <CGAL/tags.h>
#include <CGAL/Uncertain.h>
#include <CGAL/certified_numeric_predicates.h>
#include <CGAL/Quotient.h>
#include <CGAL/certified_quotient_predicates.h>
#include <CGAL/Unfiltered_predicate_adaptor.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/tuple/tuple.hpp>
#include <boost/optional/optional.hpp>

#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#  include<string>
#  include<iostream>
#  include<sstream>
#  include<iomanip>
#  define CGAL_STSKEL_TRAITS_TRACE(m) \
     { \
       std::ostringstream ss ; \
       ss << std::setprecision(19) << m << std::ends ; \
       std::string s = ss.str(); \
       Straight_skeleton_traits_external_trace(s); \
     }
#else
#  define CGAL_STSKEL_TRAITS_TRACE(m)
#endif

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i {

using boost::optional ;

// boost::make_optional is provided in Boost >= 1.34, but not before, so we define our own versions here.
template<class T> optional<T> cgal_make_optional( T const& v ) { return optional<T>(v) ; }
template<class T> optional<T> cgal_make_optional( bool cond, T const& v ) { return cond ? optional<T>(v) : optional<T>() ; }

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
// This number type is provided because unlike Quotient<> is allows you to create it
// with a zero denominator.
//
template<class NT>
class Rational
{
  public:

    Rational( NT aN, NT aD ) : mN(aN), mD(aD) {}

    NT n() const { return mN ; }
    NT d() const { return mD ; }

    CGAL::Quotient<NT> to_quotient() const { return CGAL::Quotient<NT>(mN,mD) ; }

  private:

    NT mN, mD ;
} ;


//
// A straight skeleton event is defined by 3 oriented straight line segments .
// Such a triple of segments is encapsulated in this record.
//
// NOTE: A triedge is not suffixed C2 because it holds 3 segments which can be cartesian or homogeneous
//
template<class K>
class Triedge_2
{
  public:

    typedef typename K::Segment_2 Segment_2 ;

    Triedge_2( Segment_2 const& aE0, Segment_2 const& aE1, Segment_2 const& aE2 ) : mE0(aE0), mE1(aE1), mE2(aE2) {}

    Segment_2 const& e0() const { return mE0 ; }
    Segment_2 const& e1() const { return mE1 ; }
    Segment_2 const& e2() const { return mE2 ; }

    Segment_2 const& e( int idx ) const { return idx == 0 ? mE0 : idx == 1 ? mE1 : mE2 ; }

    friend std::ostream& operator << ( std::ostream& os, Triedge_2<K> const& aTriedge )
    {
      return os << "Triedge(" << aTriedge.e0() << "\n," << aTriedge.e1() << "\n," << aTriedge.e2() << ')' ;
    }

  private:

    Segment_2 mE0, mE1, mE2 ;
} ;

//
// Most calculations need to know if the are collinear edges in a triedge.
// To that effect, a triedge is "sorted" such that collinear edges, if any, are stored in e0 and e1; with e2
// being the non-collinear edge.
// The total number of collinear edges is also recorded. It cannot be >= 3.
//
// A SortedTriedge is "indeterminate" if the collinearity of the edges couldn't be determined reliably.
//
template<class K>
class Sorted_triedge_2 : public Triedge_2<K>
{
  public:

    typedef Triedge_2<K> Base ;

    typedef typename Base::Segment_2 Segment_2 ;

    Sorted_triedge_2( Segment_2 const& aE0
                    , Segment_2 const& aE1
                    , Segment_2 const& aE2
                    , int              aCollinearCount
                   )
     : Base(aE0,aE1,aE2)
     , mCCount(aCollinearCount)
     {}

    bool is_indeterminate() const { return mCCount == -1 ; }
    
    int  collinear_count() const { return mCCount  ; }

  private:

    int mCCount ;
} ;

template<class K>
struct Functor_base_2
{
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  
  typedef Triedge_2<K>        Triedge_2 ;
  typedef Sorted_triedge_2<K> Sorted_triedge_2 ;
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

  typedef typename Source_kernel::Segment_2 Source_segment_2 ;
  typedef typename Target_kernel::Segment_2 Target_segment_2 ;

  typedef Triedge_2<Source_kernel> Source_triedge_2 ;
  typedef Triedge_2<Target_kernel> Target_triedge_2 ;

  typedef boost::optional<Source_FT>      Source_optional_FT ;
  typedef boost::optional<Source_point_2> Source_optional_point_2 ;
  
  typedef boost::optional<Target_FT>      Target_optional_FT ;
  typedef boost::optional<Target_point_2> Target_optional_point_2 ;
  
  typedef boost::tuple<Source_optional_FT,Source_optional_point_2> Source_time_and_point ;
  typedef boost::tuple<Target_optional_FT,Target_optional_point_2> Target_time_and_point ;
  
  Target_FT        cvtn(Source_FT n) const  { return this->Converter::operator()(n); }

  Target_point_2   cvtp(Source_point_2 const& p) const  { return this->Converter::operator()(p); }

  Target_segment_2 cvts( Source_segment_2 const& e) const { return Target_segment_2(cvtp(e.source()), cvtp(e.target()) ) ; }
  
  Target_optional_FT cvton( Source_optional_FT n ) const
  {
    if ( n )
         return Target_optional_FT(cvtn(*n));
    else return Target_optional_FT(); 
  }
  
  Target_optional_point_2 cvtop( Source_optional_point_2 const& p ) const
  {
    if ( p )
         return Target_optional_point_2(cvtp(*p));
    else return Target_optional_point_2(); 
  }
  
  
  Target_FT        operator()(Source_FT n) const { return cvtn(n) ; }

  Target_point_2   operator()( Source_point_2 const& p) const { return cvtp(p) ; }

  Target_segment_2 operator()( Source_segment_2 const& s) const { return cvts(s); }
  
  Target_triedge_2 operator()( Source_triedge_2 const& t) const
  {
    return Target_triedge_2(cvts(t.e0()), cvts(t.e1()), cvts(t.e2()) ) ;
  }
  
  Target_time_and_point operator() ( Source_time_and_point const& v ) const
  {
    Source_optional_FT      t ;
    Source_optional_point_2 p ;
    boost::tie(t,p) = v ;
    return Target_time_and_point(cvton(t),cvtop(p));
  }
  
};

} // namespace CGAL_SS_i


//
// This macro defines a global functor adapter which allows users to use it in the followig ways:
//
// Given a 'Functor' provided by a given 'Traits' (or Kernel):
//
//   typedef typename CGAL::Functor<Traits>::type Functor ;
//   result r = CGAL::Functor<Traits>(traits)(a,b,c);
//
#define CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(functor) \
        template<class K> \
        typename K :: functor functor ( K const& aK ) \
        { \
          return aK.get((typename K :: functor const*)0);  \
        }


CGAL_END_NAMESPACE


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H //

// EOF //
