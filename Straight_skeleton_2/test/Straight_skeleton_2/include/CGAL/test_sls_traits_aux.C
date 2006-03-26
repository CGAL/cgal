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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/include/CGAL/_test_traits.C $
// $Id: _test_traits.C 28555 2006-02-15 18:54:04Z fcacciola $
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

template<class Traits, class triedge>
bool exist_event( Traits const&  aTraits, triedge const& aTriedge )
{
  return CGAL::Do_ss_event_exist_2<Traits>(aTraits)()(aTriedge.triple());
}

template<class Traits, class triedge>
void test_exist_event( int            i
                     , Traits const&  aTraits
                     , triedge const& aTriedge
                     , bool           aExpected
                     )
{
  report(i,aExpected==exist_event(aTraits,aTriedge), boost::str(boost::format("%s") % aTriedge) );
}

template<class Traits, class triedge>
CGAL::Comparison_result compare_events(Traits const& aTraits, triedge const& aTriedgeA, triedge const& aTriedgeB )
{
  return CGAL::Compare_ss_event_times_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple());
}

template<class Traits, class triedge>
void test_compare_events( int                     i
                        , Traits const&           aTraits
                        , triedge const&          aTriedgeA
                        , triedge const&          aTriedgeB
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

template<class Traits, class Point, class triedge>
CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,Point   const& aP
                                             ,triedge const& aTriedgeA
                                             ,triedge const& aTriedgeB
                                             )
{
  return CGAL::Compare_ss_event_distance_to_seed_2<Traits>(aTraits)()(aP,aTriedgeA.triple(),aTriedgeB.triple());
}

template<class Traits, class triedge>
CGAL::Comparison_result compare_sdist_to_seed(Traits  const& aTraits
                                             ,triedge const& aTriedgeA
                                             ,triedge const& aTriedgeB
                                             ,triedge const& aTriedgeC
                                             )
{
  return CGAL::Compare_ss_event_distance_to_seed_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple(),aTriedgeC.triple());
}

template<class Traits, class Point, class triedge>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , Point  const&           aP
                               , triedge const&          aTriedgeA
                               , triedge const&          aTriedgeB
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

template<class Traits, class triedge>
void test_compare_sdist_to_seed( int                     i
                               , Traits const&           aTraits
                               , triedge const&          aTriedgeA
                               , triedge const&          aTriedgeB
                               , triedge const&          aTriedgeC
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

template<class Traits, class triedge>
bool is_inside_offset_zone( Traits const&  aTraits, triedge const& aTriedgeA, triedge const& aTriedgeB )
{
  return CGAL::Is_ss_event_inside_offset_zone_2<Traits>(aTraits)()(aTriedgeA.triple(),aTriedgeB.triple());
}

template<class Traits, class triedge>
void test_is_inside_offset_zone( int            i
                               , Traits const&  aTraits
                               , triedge const& aTriedgeA
                               , triedge const& aTriedgeB
                               , bool           aExpected
                               )
{
  if ( !exist_event(aTraits,aTriedgeA) )
    report(i,false, boost::str(boost::format("Event A doesn't exist: %s") % aTriedgeA));
  else
    report(i,aExpected==is_inside_offset_zone(aTraits,aTriedgeA,aTriedgeB), boost::str(boost::format("%s,%s") % aTriedgeA % aTriedgeB) );
}

