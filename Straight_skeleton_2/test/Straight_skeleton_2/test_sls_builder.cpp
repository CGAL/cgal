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
#include<string> 
#include<iostream>      
#include<fstream>
#include<sstream>   


//#define CGAL_STRAIGHT_SKELETON_ENABLE_TRACE 0
//#define CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
void Straight_skeleton_external_trace ( std::string m )
{
  printf("%s\n",m.c_str());
}
#endif

#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
void Straight_skeleton_traits_external_trace ( std::string m )
{
  printf("%s\n",m.c_str());
}
#endif

#include <CGAL/test_sls_builder_types.h> 

#include <CGAL/Real_timer.h>

using namespace std ;

int    sFailed  = 0 ;
bool   sDumpEPS = false ;

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
      PolygonPtr lPoly( new Polygon_2() );
      
      int v_count ;
      in >> v_count ;
      for ( int j = 0 ; j < v_count ; ++ j )
      {
        double x,y ;
        in >> x >> y ;
        lPoly->push_back( Point(x,y) ) ;
      }
      if ( lPoly->size() >= 3 )
      {
        CGAL::Orientation expected = ( i == 0 ? CGAL::COUNTERCLOCKWISE : CGAL::CLOCKWISE ) ;
        
        double area = CGAL::polygon_area_2(lPoly->begin(),lPoly->end(),K());
        
        CGAL::Orientation orientation = area > 0 ? CGAL::COUNTERCLOCKWISE : area < 0 ? CGAL::CLOCKWISE : CGAL::COLLINEAR ;
        
        if ( orientation == expected )
             rRegion->push_back(lPoly);
        else rRegion->push_back( PolygonPtr( new Polygon_2 CGAL_make_vector(lPoly->rbegin(),lPoly->rend()) ) ) ;
      }
      else cerr << "Degenerate polygon in file " << file << endl ;
    }
  }
  else cerr << "Cannot open input file " << file << endl ;

  return rRegion ;
}

void dump_eps( SlsPtr sls, std::string eps )
{
 CGAL::Bbox_2 bbox ;

  for(Vertex_const_iterator vit = sls->vertices_begin(); vit != sls->vertices_end(); ++vit)
  {
    if( vit == sls->vertices_begin() )
	        bbox =        vit->point().bbox(); 
    else	bbox = bbox + vit->point().bbox();
  }

  double scale = 1 ; //1000 / (bbox.xmax() - bbox.xmin()) ;

  if ( scale < 1 )
     scale = 1 ;

  std::ofstream dump(eps.c_str());
	 dump << "%!PS-Adobe-2.0 EPSF-2.0\n%%BoundingBox:" 
       << static_cast<int>(std::floor(scale* bbox.xmin()-1)) 
       << " " 
       << static_cast<int>(std::floor(scale* bbox.ymin()-1)) 
       << " "
       << static_cast<int>(std::ceil(scale*bbox.xmax()+1)) 
       << " " 
       << static_cast<int>(std::ceil(scale*bbox.ymax()+1))
       << std::endl;

	 dump << "%%EndComments\n"
	         "gsave\n"
	         "1.0 setlinewidth\n"
	         "/cont { 0 0 0 setrgbcolor } bind def\n"
	         "/cont_w { 0.1 setlinewidth } bind def\n"
	         "/skel { 1 0 0 setrgbcolor } bind def\n"
	         "/skel_w { 1.0 setlinewidth } bind def\n"
	         "% stroke - x1 y1 x2 y2 E\n"
	         "/E {newpath moveto lineto stroke} bind def\n" 
          << std::endl;

  for(Face_const_iterator fit = sls->faces_begin(); fit != sls->faces_end(); ++fit)
  {
	   Halfedge_const_handle h = fit->halfedge();
	   Halfedge_const_handle done;
	   done = h;
	   do
    {
	     if(h->is_bisector())
	          dump << "skel\n";
	     else dump << "cont\n";
 	  
	     dump << scale * h->vertex()->point().x() 
           << " " 
           << scale * h->vertex()->point().y()
           << " "
		         << scale * h->opposite()->vertex()->point().x()
           << " "
           << scale * h->opposite()->vertex()->point().y() 
           << " E\n";

	     h = h->next();
	   } 
    while(h != done);
	 }

	 dump << "grestore\nshowpage" << std::endl;

	 dump.close();
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
      builder.enter_contour((*bit)->begin(),(*bit)->end());
    SlsPtr sls = builder.construct_skeleton() ;
    t.stop();
    bool ok = sls  ;
    cout << file << " : " << ( ok ? "OK" : "FAILED!" ) << " (" << t.time() << " seconds)." << endl ;
    if ( ok )
    {
      (*ok_list) << file << endl ;
      if ( sDumpEPS )
       dump_eps(sls, file + std::string(".eps"));
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
  bool print_usage = true ;

  int aidx = 1 ;
  while ( aidx < argc && argv[aidx][0] == '-' )
  {
    switch(argv[aidx][1])
    {
      case 'e' : sDumpEPS = true ; break ;
      default: cerr << "Invalid option: " << argv[aidx] << endl ; break ;
    }
    ++aidx ;
  }

  std::vector<std::string> cases ;
  
  if ( aidx + 1 < argc ) 
  {
    std::string folder(argv[aidx]);
    for ( int i = aidx + 1 ; i < argc ; ++ i )
      cases.push_back( folder + std::string(argv[i]) );
  }
      
  
  if ( cases.size() > 0  ) 
  {
    print_usage = false ;
    
    std::cout << "Testing straight skeleton\n";

    failed_list = new std::ofstream("./test_sls_failed_cases.txt");
    ok_list     = new std::ofstream("./test_sls_ok_cases.txt");

    if ( !failed_list->good() )
      cerr << "Unable to open failed_cases.txt report file." << endl ;

    if ( !ok_list->good() )
      cerr << "Unable to open failed_cases.txt report file." << endl ;

    try
    {
      std::string folder(argv[aidx]);
      for ( std::vector<std:.stringint i = aidx + 1 ; i < argc ; ++ i )
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

  if ( print_usage )
  {
    cout << "USAGE: test_sls_builder <options> folder file0 file1 ... fileN" << endl
         << "  <options>: " << endl
         << "     -e  Dumps result into an .eps file." << endl
         << "     -h       Prints this usage summary." << endl ;
  }

  return sFailed == 0 ? 0 : 1 ;
}
