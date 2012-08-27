// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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

