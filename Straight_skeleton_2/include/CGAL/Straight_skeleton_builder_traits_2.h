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
#ifndef CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H
#define CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H 1

#include <CGAL/Filtered_construction.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_aux.h>
#include <CGAL/Straight_skeleton_2/Straight_skeleton_builder_traits_2_aux.h>
#include <CGAL/predicates/Straight_skeleton_pred_ftC2.h>
#include <CGAL/constructions/Straight_skeleton_cons_ftC2.h>

CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i {

template<class K>
struct Construct_ss_triedge_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;
  
  typedef typename Base::FT        FT ;
  typedef typename Base::Segment_2 Segment_2 ;
  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Triedge_2    result_type ;
  typedef Arity_tag<3> Arity ;

  Triedge_2 operator() ( Segment_2 const& aA, Segment_2 const& aB, Segment_2 const& aC ) const
  {
    return Triedge_2(aA,aB,aC);
  }
};

template<class K>
struct Construct_ss_sorted_triedge_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;
  
  typedef typename Base::FT               FT ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Triedge_2        Triedge_2 ;
  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Triedge_2    result_type ;
  typedef Arity_tag<3> Arity ;

  Sorted_triedge_2 operator() ( Triedge_2 const& aT, Triedge_collinearity aC )
  {
    return construct_sorted_triedge(aT,aC);
  }
};

template<class K>
struct Get_ss_triedge_collinearity_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2 Triedge_2 ;

  typedef Uncertain<Triedge_collinearity> result_type ;
  typedef Arity_tag<1>                    Arity ;

  Uncertain<Triedge_collinearity> operator() ( Triedge_2 const& aTriedge ) const
  {
    Uncertain<Triedge_collinearity> rResult = certified_triedge_collinearity(aTriedge);
    
    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Get_ss_triedge_collinearity_2",aTriedge);

    return rResult ;
  }
};

template<class K>
struct Do_ss_event_exist_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<1>    Arity ;

  Uncertain<bool> operator() ( Sorted_triedge_2 const& aSTriedge ) const
  {
    Uncertain<bool> rResult = exist_offset_lines_isec2(aSTriedge) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Exist_event",aSTriedge);

    return rResult ;
  }
};

template<class K>
struct Compare_ss_event_distance_to_seed_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<3>                 Arity ;

  Uncertain<Comparison_result> operator() ( Point_2          const& aP
                                          , Sorted_triedge_2 const& aL
                                          , Sorted_triedge_2 const& aR
                                          ) const
  {
    return compare_offset_lines_isec_dist_to_pointC2(cgal_make_optional(aP),aL,aR) ;
  }

  Uncertain<Comparison_result> operator() ( Sorted_triedge_2 const& aS
                                          , Sorted_triedge_2 const& aL
                                          , Sorted_triedge_2 const& aR
                                          ) const
  {
    return compare_offset_lines_isec_dist_to_pointC2(aS,aL,aR) ;
  }

};

template<class K>
struct Compare_ss_event_times_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Uncertain<Comparison_result> result_type ;
  typedef Arity_tag<2>                 Arity ;

  Uncertain<Comparison_result> operator() ( Sorted_triedge_2 const& aL, Sorted_triedge_2 const& aR ) const
  {
    Uncertain<Comparison_result> rResult = compare_offset_lines_isec_timesC2(aL,aR) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Compare_event_times","L: " << aL << "\nR:" << aR );

    return rResult ;
  }
};

template<class K>
struct Is_ss_event_inside_offset_zone_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Triedge_2        Triedge_2 ;
  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Sorted_triedge_2 const& aE, Triedge_2 const& aZ ) const
  {
    Uncertain<bool> rResult = is_offset_lines_isec_inside_offset_zoneC2(aE,aZ) ;

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Is_event_inside_offset_zone","E=" << aE << "\nZ=" << aZ);

    return rResult ;
  }
};

template<class K>
struct Are_ss_events_simultaneous_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Sorted_triedge_2 const& aA, Sorted_triedge_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_events_simultaneousC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_events_simultaneous","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_collinear_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_collinearC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_collinear","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Are_ss_edges_parallel_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::Segment_2 Segment_2 ;

  typedef Uncertain<bool> result_type ;
  typedef Arity_tag<2>    Arity ;

  Uncertain<bool> operator() ( Segment_2 const& aA, Segment_2 const& aB ) const
  {
    Uncertain<bool> rResult = are_edges_parallelC2(aA,aB);

    CGAL_STSKEL_ASSERT_PREDICATE_RESULT(rResult,K,"Are_ss_edges_parallel","A=" << aA << "\nB=" << aB);

    return rResult ;
  }
};

template<class K>
struct Construct_ss_event_time_and_point_2 : Functor_base_2<K>
{
  typedef Functor_base_2<K> Base ;

  typedef typename Base::FT               FT ;
  typedef typename Base::Point_2          Point_2 ;
  typedef typename Base::Segment_2        Segment_2 ;
  typedef typename Base::Triedge_2        Triedge_2 ;
  typedef typename Base::Sorted_triedge_2 Sorted_triedge_2 ;

  typedef boost::tuple<FT,Point_2> rtype ;
  
  typedef boost::optional<rtype> result_type ;
  
  typedef Arity_tag<1>  Arity ;

  result_type operator() ( Sorted_triedge_2 const& triedge ) const
  {
    bool lOK = false ;
    
    FT      t(0) ;
    Point_2 i = ORIGIN ;

    CGAL_assertion(triedge.collinear_count() < 3) ;

    optional< Rational<FT> > ot = compute_offset_lines_isec_timeC2(triedge);

    if ( !!ot && certainly( CGAL_NTS certified_is_not_zero(ot->d()) ) )
    {
      t = ot->n() / ot->d();
      
      optional<Point_2> oi = construct_offset_lines_isecC2(triedge);
      if ( oi )
      {
        i = *oi ;
        lOK = is_point_calculation_accurate(t,i,triedge);
      } 
    }
    
    CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(lOK,K,"Construct_ss_event_time_and_point_2",triedge);


    return cgal_make_optional(lOK,boost::make_tuple(t,i)) ;
  }
    
  bool is_point_calculation_accurate( double t, Point_2 const& p, Triedge_2 const& triedge ) const 
  {
    Segment_2 const& e0 = triedge.e0() ; 
    Segment_2 const& e1 = triedge.e1() ;
    Segment_2 const& e2 = triedge.e2() ;
    
    Point_2 const& e0s = e0.source();
    Point_2 const& e0t = e0.target();
    
    Point_2 const& e1s = e1.source();
    Point_2 const& e1t = e1.target();
    
    Point_2 const& e2s = e2.source();
    Point_2 const& e2t = e2.target();
  
    double d0 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e0s.x(),e0s.y(),e0t.x(),e0t.y());
    double d1 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e1s.x(),e1s.y(),e1t.x(),e1t.y());
    double d2 = squared_distance_from_point_to_lineC2(p.x(),p.y(),e2s.x(),e2s.y(),e2t.x(),e2t.y());
  
    double tt = CGAL_NTS square(t);
    
    double diff0 = CGAL_NTS square(d0-tt);
    double diff1 = CGAL_NTS square(d1-tt);
    double diff2 = CGAL_NTS square(d2-tt);
    
    double const zero = 1e-16 ;

    return ( diff0 <= zero && diff1 <= zero && diff2 <= zero ) ;
  }

  template<class NT>  
  bool is_point_calculation_accurate( NT const& t, Point_2 const& p, Triedge_2 const& triedge ) const 
  {
    return true ;
  }

};

} // namespace CGAL_SS_i

template<class K>
struct Straight_skeleton_builder_traits_2_functors
{
  typedef CGAL_SS_i::Get_ss_triedge_collinearity_2      <K> Get_ss_triedge_collinearity_2 ;
  typedef CGAL_SS_i::Do_ss_event_exist_2                <K> Do_ss_event_exist_2 ;
  typedef CGAL_SS_i::Compare_ss_event_times_2           <K> Compare_ss_event_times_2 ;
  typedef CGAL_SS_i::Compare_ss_event_distance_to_seed_2<K> Compare_ss_event_distance_to_seed_2 ;
  typedef CGAL_SS_i::Is_ss_event_inside_offset_zone_2   <K> Is_ss_event_inside_offset_zone_2 ;
  typedef CGAL_SS_i::Are_ss_events_simultaneous_2       <K> Are_ss_events_simultaneous_2 ;
  typedef CGAL_SS_i::Are_ss_edges_parallel_2            <K> Are_ss_edges_parallel_2 ;
  typedef CGAL_SS_i::Are_ss_edges_collinear_2           <K> Are_ss_edges_collinear_2 ;
  typedef CGAL_SS_i::Construct_ss_event_time_and_point_2<K> Construct_ss_event_time_and_point_2 ;
  typedef CGAL_SS_i::Construct_ss_triedge_2             <K> Construct_ss_triedge_2 ;
  typedef CGAL_SS_i::Construct_ss_sorted_triedge_2      <K> Construct_ss_sorted_triedge_2 ;
} ;

template<class K>
struct Straight_skeleton_builder_traits_2_base
{
  typedef K Kernel ;
  
  typedef typename K::FT        FT ;
  typedef typename K::Point_2   Point_2 ;
  typedef typename K::Segment_2 Segment_2 ;
  
  typedef CGAL_SS_i::Triedge_2       <K> Triedge_2 ;
  typedef CGAL_SS_i::Sorted_triedge_2<K> Sorted_triedge_2 ;

  template<class F> F get( F const* = 0 ) const { return F(); }
  
} ;


template<class Is_filtered_kernel, class K> class Straight_skeleton_builder_traits_2_impl ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_false,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef Straight_skeleton_builder_traits_2_functors<K> Unfiltering ;

public:

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Get_ss_triedge_collinearity_2>
    Get_ss_triedge_collinearity_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Do_ss_event_exist_2>
    Do_ss_event_exist_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_times_2>
    Compare_ss_event_times_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Compare_ss_event_distance_to_seed_2>
    Compare_ss_event_distance_to_seed_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Is_ss_event_inside_offset_zone_2>
    Is_ss_event_inside_offset_zone_2 ;

  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_events_simultaneous_2>
    Are_ss_events_simultaneous_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_parallel_2>
    Are_ss_edges_parallel_2 ;
    
  typedef Unfiltered_predicate_adaptor<typename Unfiltering::Are_ss_edges_collinear_2>
    Are_ss_edges_collinear_2 ;

  typedef typename Unfiltering::Construct_ss_event_time_and_point_2 Construct_ss_event_time_and_point_2 ;
  typedef typename Unfiltering::Construct_ss_triedge_2              Construct_ss_triedge_2 ;
  typedef typename Unfiltering::Construct_ss_sorted_triedge_2       Construct_ss_sorted_triedge_2 ;
} ;

template<class K>
class Straight_skeleton_builder_traits_2_impl<Tag_true,K> : public Straight_skeleton_builder_traits_2_base<K>
{
  typedef typename K::EK EK ;
  typedef typename K::FK FK ;

  typedef Straight_skeleton_builder_traits_2_functors<EK> Exact ;
  typedef Straight_skeleton_builder_traits_2_functors<FK> Filtering ;
  typedef Straight_skeleton_builder_traits_2_functors<K>  Unfiltering ;

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

  typedef Filtered_predicate<typename Exact    ::Get_ss_triedge_collinearity_2
                            ,typename Filtering::Get_ss_triedge_collinearity_2
                            , C2E
                            , C2F
                            >
                            Get_ss_triedge_collinearity_2 ;
                            
  typedef Filtered_predicate<typename Exact    ::Do_ss_event_exist_2
                            ,typename Filtering::Do_ss_event_exist_2
                            , C2E
                            , C2F
                            >
                            Do_ss_event_exist_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_times_2
                            , typename Filtering::Compare_ss_event_times_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_times_2 ;

  typedef Filtered_predicate< typename Exact    ::Compare_ss_event_distance_to_seed_2
                            , typename Filtering::Compare_ss_event_distance_to_seed_2
                            , C2E
                            , C2F
                            >
                            Compare_ss_event_distance_to_seed_2 ;
    
  typedef Filtered_predicate< typename Exact    ::Is_ss_event_inside_offset_zone_2
                            , typename Filtering::Is_ss_event_inside_offset_zone_2
                            , C2E
                            , C2F
                            >
                            Is_ss_event_inside_offset_zone_2 ;

  typedef Filtered_predicate< typename Exact    ::Are_ss_events_simultaneous_2
                            , typename Filtering::Are_ss_events_simultaneous_2
                            , C2E
                            , C2F
                            >
                            Are_ss_events_simultaneous_2 ;

  typedef Filtered_predicate< typename Exact    ::Are_ss_edges_parallel_2
                            , typename Filtering::Are_ss_edges_parallel_2
                            , C2E
                            , C2F
                            >
                            Are_ss_edges_parallel_2 ;
                            
  typedef Filtered_predicate< typename Exact    ::Are_ss_edges_collinear_2
                            , typename Filtering::Are_ss_edges_collinear_2
                            , C2E
                            , C2F
                            >
                            Are_ss_edges_collinear_2 ;
                            
  typedef CGAL_SS_i::Exceptionless_filtered_construction< typename Unfiltering::Construct_ss_event_time_and_point_2
                                                        , typename Exact      ::Construct_ss_event_time_and_point_2
                                                        , typename Unfiltering::Construct_ss_event_time_and_point_2 
                                                        , C2E
                                                        , C2C 
                                                        , E2C
                                                        , C2C 
                                                        >
                                                        Construct_ss_event_time_and_point_2 ;

  typedef typename Unfiltering::Construct_ss_triedge_2        Construct_ss_triedge_2 ;
  typedef typename Unfiltering::Construct_ss_sorted_triedge_2 Construct_ss_sorted_triedge_2 ;

} ;

template<class K>
class Straight_skeleton_builder_traits_2
  : public Straight_skeleton_builder_traits_2_impl<typename CGAL_SS_i::Is_filtering_kernel<K>::type, K>
{
} ;

CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Get_ss_triedge_collinearity_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Do_ss_event_exist_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_times_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Compare_ss_event_distance_to_seed_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Is_ss_event_inside_offset_zone_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_events_simultaneous_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_parallel_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Are_ss_edges_collinear_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_event_time_and_point_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_triedge_2);
CGAL_STRAIGHT_SKELETON_CREATE_FUNCTOR_ADAPTER(Construct_ss_sorted_triedge_2);

CGAL_END_NAMESPACE

#endif // CGAL_STRAIGHT_SKELETON_BUILDER_TRAITS_2_H //
// EOF //
