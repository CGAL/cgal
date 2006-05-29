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
#ifndef CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H 1

#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Polygon_offset_pred_ftC2.h>
#include <CGAL/constructions/Polygon_offset_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i {

template<class K>
struct Compare_offset_against_event_time_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT      FT ;
  typedef typename Base::Edge    Edge ;
  typedef typename Base::Triedge Triedge ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( FT aT, Triedge const& aE ) const
  {
    return compare_offset_against_isec_timeC2(aT,aE) ;
  }
};


template<class K>
struct Construct_offset_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT      FT ;
  typedef typename Base::Point_2 Point_2 ;
  typedef typename Base::Vertex  Vertex ;
  typedef typename Base::Edge    Edge ;

  typedef boost::optional< Point_2 > result_type ;
  
  typedef Arity_tag<3> Arity ;

  boost::optional<Point_2> operator() ( FT aT, Edge const& aE0, Edge const& aE1 ) const
  {
    FT ix(0.0), iy(0.0);
    
    optional<Vertex> i = construct_offset_pointC2(aT,aE0,aE1);
    
    if ( i )
    {
      ix = i->x();
      iy = i->y();
    }

    return make_optional(i, Point_2(ix,iy) ) ;
  }
};


} // namespace CGAL_SS_i

template<class K>
struct Polygon_offset_builder_traits_2_functors
{
  typedef CGAL_SS_i::Compare_offset_against_event_time_2<K> Compare_offset_against_event_time_2 ;
  typedef CGAL_SS_i::Construct_offset_point_2           <K> Construct_offset_point_2 ;
  typedef CGAL_SS_i::Construct_ss_edge_2                <K> Construct_ss_edge_2 ;
  typedef CGAL_SS_i::Construct_ss_triedge_2             <K> Construct_ss_triedge_2 ;
} ;

template<class K>
struct Polygon_offset_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT      FT ;
  typedef typename K::Point_2 Point_2 ;

  template<class F> F get( F const* = 0 ) const { return F(); }
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
  typedef typename Unfiltering::Construct_ss_edge_2      Construct_ss_edge_2 ;
  typedef typename Unfiltering::Construct_ss_triedge_2   Construct_ss_triedge_2 ;

} ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_true,K> : public Polygon_offset_builder_traits_2_base<K>
{
  typedef Polygon_offset_builder_traits_2_functors<typename K::EK> Exact ;
  typedef Polygon_offset_builder_traits_2_functors<typename K::FK> Filtering ;
  typedef Polygon_offset_builder_traits_2_functors<K>              Unfiltering ;

  typedef typename K::C2E BaseC2E ;
  typedef typename K::C2F BaseC2F ;

  typedef CGAL_SS_i::Triedge_converter<BaseC2E> C2E ;
  typedef CGAL_SS_i::Triedge_converter<BaseC2F> C2F ;

public:

  typedef Filtered_predicate<typename Exact    ::Compare_offset_against_event_time_2
                            ,typename Filtering::Compare_offset_against_event_time_2
                            , C2E
                            , C2F
                            >
                            Compare_offset_against_event_time_2 ;

  typedef typename Unfiltering::Construct_offset_point_2 Construct_offset_point_2 ;
  typedef typename Unfiltering::Construct_ss_edge_2      Construct_ss_edge_2 ;
  typedef typename Unfiltering::Construct_ss_triedge_2   Construct_ss_triedge_2 ;
} ;

template<class K>
class Polygon_offset_builder_traits_2
  : public Polygon_offset_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K >
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_offset_against_event_time_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_offset_point_2);

CGAL_END_NAMESPACE


#endif // CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H //
// EOF //
