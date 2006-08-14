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

Seeded_trisegment create_seeded_trisegment( Traits const&  aTraits, triple const& aTriple )
{
  Trisegment const& lTrisegment = aTriple.trisegment();
  Trisegment lseed, rseed ;
  return CGAL::Construct_ss_seeded_trisegment_2(aTraits)(lTrisegment,lseed,rseed);
}


bool exist_event( Traits const&  aTraits, triple const& aTriple )
{
  return CGAL::Do_ss_event_exist_2(aTraits)(create_seeded_trisegment(aTraits,aTriple));
}

template<class Traits, class triple>
void test_exist_event( int            i
                     , Traits const&  aTraits
                     , triple const& aTriple
                     , bool           aExpected
                     )
{
  report(i,aExpected==exist_event(aTraits,aTriple), boost::str(boost::format("%s") % aTriple) );
}

CGAL::Comparison_result compare_events(Traits const& aTraits, triple const& aTripleA, triple const& aTripleB )
{
  return CGAL::Compare_ss_event_times_2(aTraits)(create_seeded_trisegment(aTraits,aTripleA),create_seeded_trisegment(aTraits,aTripleB));
}

template<class Traits, class triple>
void test_compare_events( int                     i
                        , Traits const&           aTraits
                        , triple const&          aTripleA
                        , triple const&          aTripleB
                        , CGAL::Comparison_result aExpected
                        )
{
  if ( !exist_event(aTraits,aTripleA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTripleA));
  else if ( !exist_event(aTraits,aTripleB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTripleB));
  else
    report(i,aExpected==compare_events(aTraits,aTripleA,aTripleB), boost::str(boost::format("%s,%s") % aTripleA % aTripleB) );
}

CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,Point   const& aP
                                             ,triple const& aTripleA
                                             ,triple const& aTripleB
                                             )
{
  return CGAL::Compare_ss_event_distance_to_seed_2(aTraits)(aP,create_seeded_trisegment(aTraits,aTripleA),create_seeded_trisegment(aTraits,aTripleB));
}

template<class Traits, class triple>
CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,triple const& aTripleA
                                             ,triple const& aTripleB
                                             ,triple const& aTripleC
                                             )
{
  return CGAL::Compare_ss_event_distance_to_seed_2(aTraits)(create_seeded_trisegment(aTraits,aTripleA)
                                                           ,create_seeded_trisegment(aTraits,aTripleB)
                                                           ,create_seeded_trisegment(aTraits,aTripleC)
                                                           );
}

template<class Traits, class Point, class triple>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , Point  const&           aP
                               , triple const&           aTripleA
                               , triple const&           aTripleB
                               , CGAL::Comparison_result aExpected
                               )
{
  if ( !exist_event(aTraits,aTripleA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTripleA));
  else if ( !exist_event(aTraits,aTripleB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTripleB));
  else
    report(i
          ,aExpected==compare_sdist_to_seed(aTraits,aP,aTripleA,aTripleB)
          ,boost::str(boost::format("%s ; %s,%s") % aP % aTripleA % aTripleB)
          );
}

template<class Traits, class triple>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , triple const&           aTripleA
                               , triple const&           aTripleB
                               , triple const&           aTripleC
                               , CGAL::Comparison_result aExpected
                               )
{
  if ( !exist_event(aTraits,aTripleA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTripleA));
  else if ( !exist_event(aTraits,aTripleB) )
    report(i,false, boost::str(boost::format("Event B doesn't exist: %s") % aTripleB));
  else if ( !exist_event(aTraits,aTripleC) )
    report(i,false, boost::str(boost::format("Event C doesn't exist: %s") % aTripleC));
  else
    report(i
          ,aExpected==compare_sdist_to_seed(aTraits,aTripleA,aTripleB,aTripleC)
          ,boost::str(boost::format("%s,%s,%s") % aTripleA % aTripleB % aTripleC)
          );
}

bool is_inside_offset_zone( Traits const&  aTraits, triple const& aTripleA, triple const& aTripleB )
{
  return CGAL::Is_ss_event_inside_offset_zone_2(aTraits)(create_seeded_trisegment(aTraits,aTripleA),create_seeded_trisegment(aTraits,aTripleB));
}

template<class Traits, class triple>
void test_is_inside_offset_zone( int            i
                               , Traits const&  aTraits
                               , triple const& aTripleA
                               , triple const& aTripleB
                               , bool           aExpected
                               )
{
  if ( !exist_event(aTraits,aTripleA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTripleA));
  else
    report(i,aExpected==is_inside_offset_zone(aTraits,aTripleA,aTripleB), boost::str(boost::format("%s,%s") % aTripleA % aTripleB) );
}

