// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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

#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
#  include<string>
#  include<iostream>
#  include<sstream>
#  include<iomanip>
#  define CGAL_SSTRAITS_TRACE(m) \
     { \
       std::ostringstream ss ; \
       ss << std::setprecision(19) << m << std::ends ; \
       std::string s = ss.str(); \
       Straight_skeleton_traits_external_trace(s); \
     }
#else
#  define CGAL_SSTRAITS_TRACE(m)
#endif

CGAL_BEGIN_NAMESPACE

namespace CGAL_SLS_i {

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

template<class FT>
class Rational // Permits zero denominator, unlike Quotient<>
{
  public:

    Rational( FT aN, FT aD ) : mN(aN), mD(aD) {}

    FT n() const { return mN ; }
    FT d() const { return mD ; }

    Quotient<FT> to_quotient() const { return Quotient<FT>(mN,mD) ; }

  private:

    FT mN, mD ;
} ;

template<class FT>
class Vertex
{
  public:

    Vertex( FT aX, FT aY ) : mX(aX), mY(aY) {}

    FT x() const { return mX ; }
    FT y() const { return mY ; }

   friend std::ostream& operator << ( std::ostream& os, Vertex<FT> const& aV )
   {
     return os << "Vertex(" << to_double(aV.x()) << ',' << to_double(aV.y()) << ')';
   }

  private:

    FT mX, mY ;
} ;

template<class FT>
class Edge
{
  public:

    typedef Vertex<FT> Vertex ;

    Edge( Vertex const& aS, Vertex const& aT ) : mS(aS), mT(aT) {}

    Vertex const& s() const { return mS ; }
    Vertex const& t() const { return mT ; }

   friend std::ostream& operator << ( std::ostream& os, Edge<FT> const& aE )
   {
     return os << "Edge(" << aE.s() << ',' << aE.t() << ')' ;
   }

  private:

    Vertex mS, mT ;
} ;

template<class FT>
class Triedge
{
  public:

    typedef Edge<FT> Edge ;

    Triedge( Edge const& aE0, Edge const& aE1, Edge const& aE2 ) : mE0(aE0), mE1(aE1), mE2(aE2) {}

    Edge const& e0() const { return mE0 ; }
    Edge const& e1() const { return mE1 ; }
    Edge const& e2() const { return mE2 ; }

    Edge const& e( int idx ) const { return idx == 0 ? mE0 : idx == 1 ? mE1 : mE2 ; }

    friend std::ostream& operator << ( std::ostream& os, Triedge<FT> const& aTriedge )
    {
      return os << "Triedge(" << aTriedge.e0() << "\n," << aTriedge.e1() << "\n," << aTriedge.e2() << ')' ;
    }

  private:

    Edge mE0, mE1, mE2 ;
} ;

template<class FT>
class SortedTriedge : public Triedge<FT>
{
  public:

    typedef Triedge<FT> Base ;

    typedef typename Base::Edge Edge ;

    SortedTriedge( Edge const& aE0
                 , Edge const& aE1
                 , Edge const& aE2
                 , bool        aIsValid
                 , bool        aIsDegenerate
                 )
     : Base(aE0,aE1,aE2)
     , mIsValid(aIsValid)
     , mIsDegenerate(aIsDegenerate)
     {}

    bool is_valid     () const { return mIsValid ; }
    bool is_degenerate() const { return mIsDegenerate ; }

  private:

    bool mIsValid, mIsDegenerate ;
} ;

template<class FT>
class Line
{
  public:

    Line( FT aA, FT aB, FT aC ) : mA(aA), mB(aB), mC(aC) {}

    FT a() const { return mA ; }
    FT b() const { return mB ; }
    FT c() const { return mC ; }

  private:

    FT mA,mB,mC ;
} ;

template<class K>
struct Sls_functor_base_2
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  typedef Vertex       <FT> Vertex ;
  typedef Edge         <FT> Edge   ;
  typedef Triedge      <FT> Triedge ;
  typedef SortedTriedge<FT> SortedTriedge ;

  static Vertex toVertex( Point_2 const& p ) { return Vertex(p.x(),p.y()) ; }

};

template<class Converter>
struct Triedge_converter : Converter
{

  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Source_kernel::FT SFT ;
  typedef typename Target_kernel::FT TFT ;

  typedef typename Source_kernel::Point_2 Source_point_2 ;
  typedef typename Target_kernel::Point_2 Target_point_2 ;

  typedef Vertex<SFT> Source_vertex ;
  typedef Vertex<TFT> Target_vertex ;

  typedef Edge<SFT> Source_edge ;
  typedef Edge<TFT> Target_edge ;

  typedef Triedge<SFT> Source_triedge ;
  typedef Triedge<TFT> Target_triedge ;

  TFT cvtn(SFT n) const  { return Converter::operator()(n); }

  Target_vertex cvtv(Source_vertex const& v) const  { return Target_vertex( cvtn(v.x()), cvtn(v.y()) ); }

  Target_point_2 cvtp(Source_point_2 const& p) const  { return Converter::operator()(p); }

  TFT operator()(SFT n) const { return cvtn(n) ; }

  Target_point_2 operator()( Source_point_2 const& p) const { return cvtp(p) ; }

  Target_triedge  operator()( Source_triedge const& s) const
  {
    return Target_triedge( Target_edge(cvtv(s.e0().s()), cvtv(s.e0().t()) )
                         , Target_edge(cvtv(s.e1().s()), cvtv(s.e1().t()) )
                         , Target_edge(cvtv(s.e2().s()), cvtv(s.e2().t()) )
                         ) ;
  }
};

} // namespace CGAL_SLS_i


//
// This macro defines a global functor adapter which allows users to use it in the followig ways:
//
// Given a 'Functor' provided by a given 'Traits' (or Kernel):
//
//   typedef typename CGAL::Functor<Traits>::type Functor ;
//   result r = CGAL::Functor<Traits>(traits)(a,b,c);
//
#define SLS_CREATE_FUNCTOR_ADAPTER(functor) \
        template<class K> \
        class functor \
        { \
          K const& mK; \
          public: \
            typedef typename K :: functor type ; \
            functor( K const& aK ) : mK(aK) {} \
            type operator()() const { return mK.get<type>(); } \
        }


CGAL_END_NAMESPACE


#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_AUX_H //

// EOF //
