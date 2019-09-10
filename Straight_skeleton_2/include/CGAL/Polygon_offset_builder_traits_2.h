// Copyright (c) 2005-2008 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
#ifndef CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H
#define CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H 1

#include <CGAL/license/Straight_skeleton_2.h>


#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/Straight_skeleton_builder_traits_2.h>
#include <CGAL/predicates/Polygon_offset_pred_ftC2.h>
#include <CGAL/constructions/Polygon_offset_cons_ftC2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

namespace CGAL {

namespace CGAL_SS_i {

template<class K>
struct Compare_offset_against_event_time_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef Uncertain<Comparison_result> result_type ;

  Uncertain<Comparison_result> operator() ( FT const& aT, Trisegment_2_ptr const& aE ) const
  {
    return compare_offset_against_isec_timeC2(aT,aE) ;
  }
};


template<class K>
struct Construct_offset_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Trisegment_2_ptr Trisegment_2_ptr ;

  typedef boost::optional<Point_2> result_type ;
  

  result_type operator() ( FT               const& aT
                         , Segment_2        const& aE0
                         , Segment_2        const& aE1 
                         , Trisegment_2_ptr const& aNode
                         ) const
  {
    typename Has_inexact_constructions<K>::type has_inexact_constructions
    CGAL_SUNPRO_INITIALIZE( = typename Has_inexact_constructions<K>::type()) ;
    
    return calc(aT, aE0, aE1, aNode, has_inexact_constructions);
  }
  
  result_type calc ( FT               const& aT
                   , Segment_2        const& aE0
                   , Segment_2        const& aE1 
                   , Trisegment_2_ptr const& aNode
                   , Tag_false        // kernel already has exact constructions
                   ) const
  {
    result_type p = construct_offset_pointC2(aT,aE0,aE1,aNode);
    
    CGAL_stskel_intrinsic_test_assertion(!p || (p && !is_point_calculation_clearly_wrong(aT,*p,aE0,aE1)));
    
    return p ;
  }
  
  result_type calc ( FT               const& aT
                   , Segment_2        const& aE0
                   , Segment_2        const& aE1 
                   , Trisegment_2_ptr const& aNode
                   , Tag_true         // kernel does not provides exact constructions
                   ) const
  {
    typedef Exact_predicates_exact_constructions_kernel EK ;
    
    typedef Cartesian_converter<K,EK> BaseC2E;
    typedef Cartesian_converter<EK,K> BaseE2C;
    
    SS_converter<BaseC2E> C2E ;  
    SS_converter<BaseE2C> E2C ;  
    
    result_type p = E2C(construct_offset_pointC2(C2E(aT),C2E(aE0),C2E(aE1),C2E(aNode)));
    
    CGAL_stskel_intrinsic_test_assertion(!p || (p && !is_point_calculation_clearly_wrong(aT,*p,aE0,aE1)));
    
    return p ;
  }
  
  bool is_point_calculation_clearly_wrong( FT const& t, Point_2 const& p, Segment_2 const& aE0, Segment_2 const& aE1 ) const
  { 
    bool rR = false ;
    
    if ( is_possibly_inexact_time_clearly_not_zero(t) )
    {
      Point_2 const& e0s = aE0.source();
      Point_2 const& e0t = aE0.target();
      
      Point_2 const& e1s = aE1.source();
      Point_2 const& e1t = aE1.target();
    
      FT const very_short(0.1);
      FT const very_short_squared = CGAL_NTS square(very_short);
      
      FT l0 = squared_distance(e0s,e0t) ;
      FT l1 = squared_distance(e1s,e1t) ;
      
      bool e0_is_not_very_short = l0 > very_short_squared ;
      bool e1_is_not_very_short = l1 > very_short_squared ;
      
      FT d0 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e0s.x(),e0s.y(),e0t.x(),e0t.y()).to_nt();
      FT d1 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e1s.x(),e1s.y(),e1t.x(),e1t.y()).to_nt();
    
      FT tt = CGAL_NTS square(t) ;
      
      bool e0_is_clearly_wrong = e0_is_not_very_short && is_possibly_inexact_distance_clearly_not_equal_to(d0,tt) ;
      bool e1_is_clearly_wrong = e1_is_not_very_short && is_possibly_inexact_distance_clearly_not_equal_to(d1,tt) ;
              
      rR = e0_is_clearly_wrong || e1_is_clearly_wrong ;        
      
      CGAL_stskel_intrinsic_test_trace_if(rR
                                        , "\nOffset point calculation is clearly wrong:"
                                          << "\ntime=" << t << " p=" << p2str(p) << " e0=" << s2str(aE0) << " e1=" << s2str(aE1)
                                          << "\nl0=" << inexact_sqrt(l0) << " l1=" << inexact_sqrt(l1)
                                          << "\nd0=" << d0 << " d1=" << d1 << " tt=" << tt 
                                        ) ;
    }
    
    return rR ;
  }
};


} // namespace CGAL_SS_i

template<class K>
struct Polygon_offset_builder_traits_2_functors
{
  typedef CGAL_SS_i::Compare_offset_against_event_time_2<K> Compare_offset_against_event_time_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2           <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Construct_offset_point_2           <K> Construct_offset_point_2 ;
  typedef CGAL_SS_i::Construct_ss_trisegment_2          <K> Construct_ss_trisegment_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2<K> Construct_ss_event_time_and_point_2 ;
} ;

template<class K>
struct Polygon_offset_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  
  typedef CGAL_SS_i::Trisegment_2<K> Trisegment_2 ;
  
  typedef typename Trisegment_2::Self_ptr Trisegment_2_ptr ;

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

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;
    
  typedef typename Unfiltering::Construct_offset_point_2            Construct_offset_point_2 ;
  typedef typename Unfiltering::Construct_ss_trisegment_2           Construct_ss_trisegment_2 ;
  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;

} ;

template<class K>
class Polygon_offset_builder_traits_2_impl<Tag_true,K> : public Polygon_offset_builder_traits_2_base<K>
{
  typedef typename K::Exact_kernel       EK ;
  typedef typename K::Approximate_kernel FK ;
  
  typedef Polygon_offset_builder_traits_2_functors<EK> Exact ;
  typedef Polygon_offset_builder_traits_2_functors<FK> Filtering ;
  typedef Polygon_offset_builder_traits_2_functors<K>  Unfiltering ;

  typedef Cartesian_converter<K,EK> BaseC2E;
  typedef Cartesian_converter<K,FK> BaseC2F;
  typedef Cartesian_converter<EK,K> BaseE2C;
  typedef Cartesian_converter<FK,K> BaseF2C;
  typedef Cartesian_converter<K,K>  BaseC2C;
  
  typedef CGAL_SS_i::SS_converter<BaseC2E> C2E ;
  typedef CGAL_SS_i::SS_converter<BaseC2F> C2F ;
  typedef CGAL_SS_i::SS_converter<BaseE2C> E2C ;
  typedef CGAL_SS_i::SS_converter<BaseF2C> F2C ;
  typedef CGAL_SS_i::SS_converter<BaseC2C> C2C ;

public:

  typedef Filtered_predicate<typename Exact    ::Compare_offset_against_event_time_2
                            ,typename Filtering::Compare_offset_against_event_time_2
                            , C2E
                            , C2F
                            >
                            Compare_offset_against_event_time_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_times_2
                            , typename Filtering::Compare_ss_event_times_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_times_2 ;
                            
  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_offset_point_2
                                                        , typename Exact      ::Construct_offset_point_2
                                                        , typename Unfiltering::Construct_offset_point_2
                                                        , C2E
                                                        , C2C
                                                        , E2C
                                                        , C2C
                                                        >
                                                        Construct_offset_point_2 ;
                                             
  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_trisegment_2
                                                        , typename Exact      ::Construct_ss_trisegment_2
                                                        , typename Unfiltering::Construct_ss_trisegment_2
                                                        , C2E
                                                        , C2C
                                                        , E2C
                                                        , C2C
                                                        >
                                                        Construct_ss_trisegment_2 ;
                                                        
  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_event_time_and_point_2
                                                        , typename Exact      ::Construct_ss_event_time_and_point_2
                                                        , typename Unfiltering::Construct_ss_event_time_and_point_2 
                                                        , C2E
                                                        , C2C 
                                                        , E2C
                                                        , C2C 
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;
} ;

template<class K>
class Polygon_offset_builder_traits_2
  : public Polygon_offset_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type,K>
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_offset_against_event_time_2)
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_offset_point_2)

} // end namespace CGAL


#endif // CGAL_POLYGON_OFFSET_BUILDER_TRAITS_2_H //
// EOF //
