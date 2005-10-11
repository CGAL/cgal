// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/_test_traits.C
// revision      : 
// revision_date : 
// author(s)     : Fernando Cacciola (fernando.cacciola@gmail.com)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================

#include<boost/format.hpp>
#include<string>

int sSucceeded = 0 ;
int sFailed    = 0 ;
bool sReportOK      = false ;
bool sReportFailed  = true ;
bool sReportSummary = true ;
    
void report( int idx, bool ok, std::string const& info = std::string("") )
{
  if (ok)
  {
    if ( sReportOK )
      std::cout << sPrefix << "Test case " << idx << " OK " << info << std::endl ;
    ++ sSucceeded ;
   }
  else 
  {
    if ( sReportFailed )
      std::cout << sPrefix << "Test case " << idx << " FAILED! " << info << std::endl ;
    ++ sFailed ;
  }
}

template<class Traits, class Triedge>
bool exist_event( Traits const&  aTraits, Triedge const& aTriedge )
{
  return CGAL::Exist_sls_event_2<Traits>(aTraits)()(aTriedge.triple());
}                    

template<class Traits, class Triedge>
void test_exist_event( int            i 
                     , Traits const&  aTraits
                     , Triedge const& aTriedge
                     , bool           aExpected 
                     )
{
  report(i,aExpected==exist_event(aTraits,aTriedge), boost::str(boost::format("%s") % aTriedge) );
}		     

template<class Traits, class Triedge>
CGAL::Comparison_result compare_events(Traits const& aTraits, Triedge const& aTriedgeA, Triedge const& aTriedgeB )
{
  return CGAL::Compare_sls_event_times_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple());
}                    

template<class Traits, class Triedge>
void test_compare_events( int                     i
                        , Traits const&           aTraits
                        , Triedge const&          aTriedgeA
                        , Triedge const&          aTriedgeB
                        , CGAL::Comparison_result aExpected 
                        )
{
  if ( !exist_event(aTraits,aTriedgeA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTriedgeA));
  else if ( !exist_event(aTraits,aTriedgeB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTriedgeB));
  else
    report(i,aExpected==compare_events(aTraits,aTriedgeA,aTriedgeB), boost::str(boost::format("%s,%s") % aTriedgeA % aTriedgeB) );
}                    

template<class Traits, class Point, class Triedge>
CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,Point   const& aP
                                             ,Triedge const& aTriedgeA
                                             ,Triedge const& aTriedgeB 
                                             )
{
  return CGAL::Compare_sls_event_distance_to_seed_2<Traits>(aTraits)()(aP,aTriedgeA.triple(),aTriedgeB.triple());
}                    

template<class Traits, class Triedge>
CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,Triedge const& aTriedgeA
                                             ,Triedge const& aTriedgeB
                                             ,Triedge const& aTriedgeC 
                                             )
{
  return CGAL::Compare_sls_event_distance_to_seed_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple(),aTriedgeC.triple());
}                    

template<class Traits, class Point, class Triedge>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , Point  const&           aP
                               , Triedge const&          aTriedgeA
                               , Triedge const&          aTriedgeB
                               , CGAL::Comparison_result aExpected 
                               )
{
  if ( !exist_event(aTraits,aTriedgeA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTriedgeA));
  else if ( !exist_event(aTraits,aTriedgeB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTriedgeB));
  else
    report(i
          ,aExpected==compare_sdist_to_seed(aTraits,aP,aTriedgeA,aTriedgeB)
          ,boost::str(boost::format("%s ; %s,%s") % aP % aTriedgeA % aTriedgeB)
          );
}                    

template<class Traits, class Triedge>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , Triedge const&          aTriedgeA
                               , Triedge const&          aTriedgeB
                               , Triedge const&          aTriedgeC
                               , CGAL::Comparison_result aExpected 
                               )
{
  if ( !exist_event(aTraits,aTriedgeA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTriedgeA));
  else if ( !exist_event(aTraits,aTriedgeB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTriedgeB));
  else if ( !exist_event(aTraits,aTriedgeC) )
    report(i,false, boost::str(boost::format("Event C doesn't exist: %s") % aTriedgeC));
  else
    report(i
          ,aExpected==compare_sdist_to_seed(aTraits,aTriedgeA,aTriedgeB,aTriedgeC)
          ,boost::str(boost::format("%s,%s,%s") % aTriedgeA % aTriedgeB % aTriedgeC)
          );
}                    
 
template<class Traits, class Triedge>
bool is_inside_offset_zone( Traits const&  aTraits, Triedge const& aTriedgeA, Triedge const& aTriedgeB )
{
  return CGAL::Is_sls_event_inside_offset_zone_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple());
}                    

template<class Traits, class Triedge>
void test_is_inside_offset_zone( int            i 
                               , Traits const&  aTraits
                               , Triedge const& aTriedgeA
                               , Triedge const& aTriedgeB
                               , bool           aExpected 
                               )
{
  if ( !exist_event(aTraits,aTriedgeA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTriedgeA));
  else
    report(i,aExpected==is_inside_offset_zone(aTraits,aTriedgeA,aTriedgeB), boost::str(boost::format("%s,%s") % aTriedgeA % aTriedgeB) );
}                    
