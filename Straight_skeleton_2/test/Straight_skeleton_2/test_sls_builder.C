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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/test_sls_traits.C $
// $Id: test_sls_traits.C 28555 2006-02-15 18:54:04Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>

#include <CGAL/test_sls_builder_types.h>

#include <CGAL/Real_timer.h>

using namespace std ;

int    sFailed = 0 ;
double sScale  = 1.0;
ofstream* failed_list = 0 ;
ofstream* ok_list     = 0 ;

RegionPtr load_region( string file )
{
  RegionPtr rRegion ;

  ifstream in(file.c_str());
  if ( in )
  {
    CGAL::set_ascii_mode(in);

    rRegion = RegionPtr( new Region() ) ;

    int ccb_count ;
    in >> ccb_count ;

    for ( int i = 0 ; i < ccb_count ; ++ i )
    {
      PolygonPtr lPoly( new Polygon() );
      in >> *lPoly;
      if ( lPoly->is_simple() )
      {
        CGAL::Orientation expected = ( i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
        if ( lPoly->orientation() != expected )
          lPoly->reverse_orientation();
        if ( sScale != 1.0 )
          lPoly = PolygonPtr( new Polygon( CGAL::transform(Transformation(CGAL::SCALING,sScale),*lPoly) ) ) ;
        rRegion->push_back(lPoly);
      }
      else cerr << "INPUT ERROR: Non-simple contour found in " << file << endl ;
    }
  }
  else cerr << "Cannot open input file " << file << endl ;

  return rRegion ;
}

void test( std::string file )
{
  RegionPtr lRegion = load_region(file);
  if ( lRegion )
  {
    CGAL::Real_timer t ;
    t.start();
    SlsBuilder builder ;
    for( Region::const_iterator bit = lRegion->begin(), ebit = lRegion->end() ; bit != ebit ; ++ bit )
      builder.enter_contour((*bit)->vertices_begin(),(*bit)->vertices_end());
    SlsPtr sls = builder.construct_skeleton() ;
    t.stop();
    bool ok = sls  ;
    cout << file << " : " << ( ok ? "OK" : "FAILED!" ) << " (" << t.time() << " seconds)." << endl ;
    if ( ok )
    {
      (*ok_list) << file << endl ;
    }
    else
    {
      (*failed_list) << file << endl ;
      ++ sFailed ;
    }
  }
}


int main( int argc, char const* argv[] )
{
  bool print_usage = false ;

  int aidx = 1 ;
  while ( aidx < argc && argv[aidx][0] == '-' )
  {
    switch(argv[aidx][1])
    {
      case 's' : sScale = atof(&argv[aidx][2]); break ;
      case 'h' : print_usage = true; break ;
      default: cerr << "Invalid option: " << argv[aidx] << endl ; break ;
    }
    ++aidx ;
  }

  if ( aidx + 1 < argc )
  {
    std::cout << "Testing Straight_skeleton_builder_2\n";

    failed_list = new std::ofstream("./sls_builder_failed_cases.txt");
    ok_list     = new std::ofstream("./sls_builder_ok_cases.txt");

    if ( !failed_list->good() )
      cerr << "Unable to open failed_cases.txt report file." << endl ;

    if ( !ok_list->good() )
      cerr << "Unable to open failed_cases.txt report file." << endl ;

    try
    {
      std::string folder(argv[aidx]);
      for ( int i = aidx + 1 ; i < argc ; ++ i )
        test(folder + std::string("/") + std::string(argv[i]) );
    }
    catch( exception x )
    {
      cerr << "Exception caught: " << x.what() << endl ;
      ++ sFailed ;
    }

    delete failed_list ;
    delete ok_list ;
  }
  else print_usage = true ;

  if ( print_usage )
  {
    cout << "USAGE: test_sls_builder <options> folder file0 file1 ... fileN" << endl
         << "  <options>: " << endl
         << "     -sSCALE  Scales each input polygon." << endl
         << "     -h       Prints this usage summary." << endl ;
  }

  return sFailed == 0 ? 0 : 1 ;
}
