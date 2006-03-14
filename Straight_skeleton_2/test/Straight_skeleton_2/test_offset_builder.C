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

#include <CGAL/test_offset_builder_types.h>

#include <CGAL/Real_timer.h>

using namespace std ;

int    sFailed = 0 ;
ofstream* failed_list = 0 ;
ofstream* ok_list     = 0 ;

PolygonPtr load_polygon( string file )
{
  PolygonPtr rPoly ;

  ifstream in(file.c_str());
  if ( in )
  {
    CGAL::set_ascii_mode(in);

    int ccb_count ;
    in >> ccb_count ; // Unused. Only the outer contour is used in this test.

    rPoly = PolygonPtr( new Polygon() );
    in >> *rPoly;
    
    if ( rPoly->orientation() != CGAL::CLOCKWISE ) // For this test we need the contour to be placed as a hole
      rPoly->reverse_orientation();
  }
  else cerr << "Cannot open input file " << file << endl ;

  return rPoly ;
}

using namespace std ;

void test( std::string file )
{
  PolygonPtr lPoly = load_polygon(file);
  if ( lPoly )
  {
     CGAL::Bbox_2 lBbox = lPoly->bbox();
     
     double w = lBbox.xmax() - lBbox.xmin();
     double h = lBbox.ymax() - lBbox.ymin();
     double s = std::sqrt(w*w+h*h);
     double lOffset = s * 0.3 ;
     
     double lMargin = compute_outer_frame_margin(lPoly->vertices_begin()
                                                ,lPoly->vertices_end  ()
                                                ,lOffset
                                                );
     double flx = lBbox.xmin() - lMargin ;
     double fhx = lBbox.xmax() + lMargin ;
     double fly = lBbox.ymin() - lMargin ;
     double fhy = lBbox.ymax() + lMargin ;
     
     Point lFrame[4]= { Point(flx,fly)
                      , Point(fhx,fly)
                      , Point(fhx,fhy)
                      , Point(flx,fhy)
                      } ;
                      
    CGAL::Real_timer t ;
    t.start();
    SlsBuilder builder ;
    builder.enter_contour(lFrame,lFrame+4);
    builder.enter_contour(lPoly->vertices_begin(),lPoly->vertices_end());
    Sls sls = builder.construct_skeleton() ;
    t.stop();
    
    bool ssok   = Sls_const_decorator(sls).is_valid(false,3);
    bool allok = false ;
    
    if ( ssok )
    {
      RegionPtr lContours( new Region ) ;
      OffsetBuilder lOffsetBuilder(sls);
      lOffsetBuilder.construct_offset_contours(lOffset, std::back_inserter(*lContours) );
      
      // Verify there are at least 2 offset contours generated.
      if ( lContours->size() > 1 )
      {
        // Find the outmost offset contour (as the one with the biggest area)
        PolygonPtr lOutmost = lContours->front();
        double lBestArea = CGAL_NTS abs (lOutmost->area());
        for( Region::const_iterator cit = successor(lContours->begin()), ecit = lContours->end() ; cit != ecit ; ++ cit )
        {
          PolygonPtr lContour = *cit ;
          if ( CGAL_NTS abs (lContour->area()) > lBestArea )
          {
            lBestArea = lContour->area();
            lOutmost  = lContour ;
          }
        }

        // Verify that the outmost offset contour is a parallelogram        
        if  ( lOutmost->size() == 5 )
        {
          double xmin = lOutmost->left_vertex  ()->x() ; 
          double xmax = lOutmost->right_vertex ()->x() ; 
          double ymin = lOutmost->bottom_vertex()->y() ; 
          double ymax = lOutmost->top_vertex   ()->y() ; 
          
          CGAL::Bbox_2 lBBox = lOutmost->bbox();
          
          // Verify that the outmost offset contour is an iso-rectangle.
          // If it is, assume this offset contour corresponds to the frame.
          if (  xmin == lBBox.xmin()
             && xmax == lBBox.xmax()
             && ymin == lBBox.ymin()
             && ymax == lBBox.ymax()
             )
          {
            double distl = xmin - flx  ;
            double distr = fhx  - xmax ;
            double distb = ymin - fly  ;
            double distt = fhy  - ymax ;
            
            double eps = 1e-5 ;
             
            // Verify the offset frame is at the right distance.
            if (  std::abs(distl-lOffset)<eps
               && std::abs(distr-lOffset)<eps
               && std::abs(distb-lOffset)<eps
               && std::abs(distt-lOffset)<eps
               )  
              allok = true ;   
          }
        }
      }
    }
        
    cout << file << " : " << ( allok ? "OK" : "FAILED!" ) << " (" << t.time() << " seconds)." << endl ;
    if ( allok )
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
      case 'h' : print_usage = true; break ;
      default: cerr << "Invalid option: " << argv[aidx] << endl ; break ;
    }
    ++aidx ;
  }

  if ( aidx + 1 < argc )
  {
    std::cout << "Testing Straight_skeleton_builder_2\n";

    failed_list = new std::ofstream("./offset_builder_failed_cases.txt");
    ok_list     = new std::ofstream("./offset_builder_ok_cases.txt");

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
    cout << "USAGE: test_offset_builder <options> folder file0 file1 ... fileN" << endl
         << "  <options>: " << endl
         << "     -h       Prints this usage summary." << endl ;
  }

  return sFailed == 0 ? 0 : 1 ;
}
