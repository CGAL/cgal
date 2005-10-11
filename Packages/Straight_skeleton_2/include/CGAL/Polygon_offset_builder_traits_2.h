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
// file          : include/CGAL/Polygon_offset_traits_2.h
// package       : Polygon_offset_2 (1.1.0)
//
// author(s)     : Fernando Cacciola
// maintainer    : Fernando Cacciola <fernando_cacciola@hotmail>
// coordinator   : Fernando Cacciola <fernando_cacciola@hotmail>
//
// ============================================================================
#ifndef CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H 1

#include <CGAL/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Polygon_offset_pred_ftC2.h>
#include <CGAL/constructions/Polygon_offset_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template<class K>
struct Compare_offset_against_event_time_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::FT          FT ;
  typedef typename Base::Edge_triple Edge_triple ;
  typedef typename Base::LineC2      LineC2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( FT aT, Edge_triple const& aE ) const
  {
    LineC2 l0,l1,l2;
    tie(l0,l1,l2) = toLineC2_triple(aE);

    return compare_offset_against_isec_timeC2(aT,l0,l1,l2) ;
  }
};


template<class K>
struct Construct_offset_point_2 : Sls_functor_base_2<K>
{
  typedef Sls_functor_base_2<K> Base ;

  typedef typename Base::FT        FT ;
  typedef typename Base::Point_2   Point_2 ;
  typedef typename Base::Edge      Edge ;
  typedef typename Base::LineC2    LineC2 ;

  typedef Point_2      result_type ;
  typedef Arity_tag<3> Arity ;

  Point_2 operator() ( FT aT, Edge const& aE0, Edge const& aE1 ) const
  {
    LineC2 l0 = toLineC2(aE0);
    LineC2 l1 = toLineC2(aE1);

    FT x,y ;
    tie(x,y) = construct_offset_pointC2(aT,l0,l1);

    return Point_2(x,y) ;
  }
};


} // namespace CGALi

template<class K>
struct Polygon_offset_builder_traits_2_functors
{
  typedef CGALi::Compare_offset_against_event_time_2<K> Compare_offset_against_event_time_2 ;
  typedef CGALi::Construct_offset_point_2           <K> Construct_offset_point_2 ;
} ;

template<class K>
struct Polygon_offset_builder_traits_2_base
{
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  template<class F> F get() const { return F(); }
} ;

template<class Is_filtered_kernel, class K> class Polygon_offset_builder_traits_2_impl ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_false,K> : public Polygon_offset_builder_traits_2_base<K>
{
  typedef Polygon_offset_builder_traits_2_functors<K> Unfiltering ;

public:

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_offset_against_event_time_2>
    Compare_offset_against_event_time_2 ;

  typedef typename Unfiltering::Construct_offset_point_2 Construct_offset_point_2 ;

} ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_true,K> : public Polygon_offset_builder_traits_2_base<K>
{
  typedef Polygon_offset_builder_traits_2_functors<typename K::EK> Exact ;
  typedef Polygon_offset_builder_traits_2_functors<typename K::FK> Filtering ;
  typedef Polygon_offset_builder_traits_2_functors<K>              Unfiltering ;

  typedef typename K::C2E BaseC2E ;
  typedef typename K::C2F BaseC2F ;

  typedef CGALi::Edge_triple_converter_2<BaseC2E> C2E ;
  typedef CGALi::Edge_triple_converter_2<BaseC2F> C2F ;

public:

  typedef Filtered_predicate<typename Exact    ::Compare_offset_against_event_time_2
                            ,typename Filtering::Compare_offset_against_event_time_2
                            , C2E
                            , C2F
                            >
                            Compare_offset_against_event_time_2 ;

  typedef typename Unfiltering::Construct_offset_point_2 Construct_offset_point_2 ;

  template<class F> F get() const { return F(); }
} ;

template<class K>
class Polygon_offset_builder_traits_2
  : public Polygon_offset_builder_traits_2_impl<typename CGALi::Is_filtering_kernel<K>::type, K >
{
} ;

SLS_CREATE_FUNCTOR_ADAPTER(Compare_offset_against_event_time_2);
SLS_CREATE_FUNCTOR_ADAPTER(Construct_offset_point_2);

CGAL_END_NAMESPACE


#endif // CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H //
// EOF //
