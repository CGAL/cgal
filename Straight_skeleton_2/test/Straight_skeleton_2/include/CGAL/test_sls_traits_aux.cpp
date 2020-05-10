// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

bool exist_event( Traits const&  aTraits, triple const& aTriple )
{
  boost::optional<FT> lMaxTime ;
  return CGAL::Do_ss_event_exist_2(aTraits)(aTriple.trisegment(), lMaxTime );
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
  return CGAL::Compare_ss_event_times_2(aTraits)(aTripleA.trisegment(),aTripleB.trisegment());
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

