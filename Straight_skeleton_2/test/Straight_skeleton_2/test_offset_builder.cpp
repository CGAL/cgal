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

    rPoly = PolygonPtr( new Polygon_2() );

    int v_count ;
    in >> v_count ;
    for ( int j = 0 ; j < v_count ; ++ j )
    {
      double x,y ;
      in >> x >> y ;
      rPoly->push_back( Point(x,y) ) ;
    }
    if ( rPoly->size() >= 3 )
    {
      double area = CGAL::polygon_area_2(rPoly->begin(),rPoly->end(),K());
        
      CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR ;
        
      if ( orientation != CGAL::CLOCKWISE )
        rPoly = PolygonPtr( new Polygon_2 CGAL_make_vector(rPoly->rbegin(),rPoly->rend()) )  ;
    }
    else
      rPoly = PolygonPtr();     
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
     bool allok = false ;
     
     CGAL::Real_timer t ;
     
     CGAL::Bbox_2 lBbox = CGAL::bbox_2(lPoly->begin(),lPoly->end());
     
     double w = lBbox.xmax() - lBbox.xmin();
     double h = lBbox.ymax() - lBbox.ymin();
     double s = std::sqrt(w*w+h*h);
     double lOffset = s * 0.3 ;
     
     boost::optional<double> lMargin = CGAL::compute_outer_frame_margin(lPoly->begin()
                                                                       ,lPoly->end  ()
                                                                       ,lOffset
                                                                       );
     if ( lMargin )
     {
       double flx = lBbox.xmin() - *lMargin ;
       double fhx = lBbox.xmax() + *lMargin ;
       double fly = lBbox.ymin() - *lMargin ;
       double fhy = lBbox.ymax() + *lMargin ;
       
       Point lFrame[4]= { Point(flx,fly)
                        , Point(fhx,fly)
                        , Point(fhx,fhy)
                        , Point(flx,fhy)
                        } ;
                        
       t.start();
       SlsBuilder builder ;
       builder.enter_contour(lFrame,lFrame+4);
       builder.enter_contour(lPoly->begin(),lPoly->end());
       SlsPtr sls = builder.construct_skeleton() ;
       t.stop();
       
       if ( sls )
       {
         RegionPtr lContours( new Region ) ;
         OffsetBuilder lOffsetBuilder(*sls);
         lOffsetBuilder.construct_offset_contours(lOffset, std::back_inserter(*lContours) );
         
         // Verify there are at least 2 offset contours generated.
         if ( lContours->size() > 1 )
         {
           // Find the outmost offset contour (as the one with the biggest area)
           PolygonPtr lOutmost = lContours->front();
           double lBestArea = CGAL_NTS abs ( CGAL::polygon_area_2(lOutmost->begin(),lOutmost->end(),K()));
           for( Region::const_iterator cit = CGAL::successor(lContours->begin()), ecit = lContours->end() ; cit != ecit ; ++ cit )
           {
             PolygonPtr lContour = *cit ;
             double lArea = CGAL_NTS abs ( CGAL::polygon_area_2(lContour->begin(),lContour->end(),K()) ) ;
             if ( lArea > lBestArea )
             {
               lBestArea = lArea ;
               lOutmost  = lContour ;
             }
           }
   
           // Verify that the outmost offset contour is a parallelogram        
           if  ( lOutmost->size() == 5 )
           {
             double xmin = CGAL::left_vertex_2  (lOutmost->begin(),lOutmost->end())->x() ; 
             double xmax = CGAL::right_vertex_2 (lOutmost->begin(),lOutmost->end())->x() ; 
             double ymin = CGAL::bottom_vertex_2(lOutmost->begin(),lOutmost->end())->y() ; 
             double ymax = CGAL::top_vertex_2   (lOutmost->begin(),lOutmost->end())->y() ; 
             
             CGAL::Bbox_2 lBBox = CGAL::bbox_2(lOutmost->begin(),lOutmost->end());
             
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
        test(folder + std::string(argv[i]) );
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
