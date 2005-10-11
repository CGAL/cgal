// ============================================================================
//
// Copyright (c) 1997-2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Straight_skeleton_builder_traits_aux_2.h
// package       : Straight_skeleton_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
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

// These are ubiquotous and will be part of the std so we use them unqualified
using boost::tuple ;
using boost::tie ;
using boost::make_tuple ;

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

namespace CGALi {

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

template<class Converter>
struct Edge_triple_converter_2 : Converter
{
  typedef typename Converter::Source_kernel Source_kernel;
  typedef typename Converter::Target_kernel Target_kernel;

  typedef typename Source_kernel::FT SFT ;
  typedef typename Target_kernel::FT TFT ;

  typedef typename Source_kernel::Point_2 Source_point_2 ;
  typedef typename Target_kernel::Point_2 Target_point_2 ;

  typedef tuple<Source_point_2, Source_point_2> Source_edge ;
  typedef tuple<Target_point_2, Target_point_2> Target_edge ;

  typedef tuple<Source_edge,Source_edge,Source_edge> Source_edge_triple ;
  typedef tuple<Target_edge,Target_edge,Target_edge> Target_edge_triple ;

  TFT cvtn(SFT n) const  { return Converter::operator()(n); }

  Target_point_2 cvtp(Source_point_2 const& p) const  { return Converter::operator()(p); }

  TFT operator()(SFT n) const { return cvtn(n) ; }

  Target_point_2 operator()( Source_point_2 const& p) const { return cvtp(p) ; }

  Target_edge operator()( Source_edge const& e) const
  {
    Source_point_2 es, et ;
    tie(es,et) = e ;

    return make_tuple(cvtp(es), cvtp(et));
  }

  Target_edge_triple  operator()( Source_edge_triple const& s) const
  {
    Source_edge s0,s1,s2;
    tie(s0,s1,s2) = s ;

    Source_point_2 s0s, s0t, s1s, s1t, s2s, s2t ;
    tie(s0s,s0t) = s0 ;
    tie(s1s,s1t) = s1 ;
    tie(s2s,s2t) = s2 ;

    return make_tuple( make_tuple(cvtp(s0s), cvtp(s0t))
                     , make_tuple(cvtp(s1s), cvtp(s1t))
                     , make_tuple(cvtp(s2s), cvtp(s2t))
                     ) ;
  }
};

template<class K>
struct Sls_functor_base_2
{
  protected :

  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  typedef tuple<FT,FT>    VertexC2 ;
  typedef tuple<FT,FT,FT> LineC2   ;

  typedef tuple<Point_2,Point_2> Edge ;

  typedef tuple<Edge,Edge,Edge>       Edge_triple ;
  typedef tuple<LineC2,LineC2,LineC2> LineC2_triple  ;

  static VertexC2 toVertexC2( Point_2 const& p ) { return make_tuple(p.x(),p.y()) ; }

  static LineC2 toLineC2( Edge const& aE )
  {
    Point_2 s,t ;
    tie(s,t) = aE ;
    return compute_normalized_line_ceoffC2(toVertexC2(s),toVertexC2(t));
  }

  static LineC2_triple toLineC2_triple( Edge_triple const& t )
  {
    Edge e0,e1,e2 ;
    tie(e0,e1,e2) = t ;
    return make_tuple(toLineC2(e0), toLineC2(e1), toLineC2(e2) ) ;
  }

  friend std::ostream& operator << ( std::ostream& os, Edge const& aEdge )
  {
    Point_2 sp, ep ; tie(sp,ep) = aEdge ;
    return os << '(' << sp.x() << ',' << sp.y() << ")->(" << ep.x() << ',' << ep.y() << ')';
  }

  friend std::ostream& operator << ( std::ostream& os, Edge_triple const& aTriple )
  {
    Edge e0,e1,e2 ; tie(e0,e1,e2) = aTriple ;
    return os << "\ne0:" << e0 << "\ne1:" << e1 << "\ne2:" << e2 ;
  }
};

} // namespace CGALi


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
