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
      PolygonPtr lPoly( new Polygon() );
      
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
        else rRegion->push_back( PolygonPtr( new Polygon(lPoly->rbegin(),lPoly->rend()) ) ) ;
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
       << scale* bbox.xmin()-1 
       << " " 
       << scale* bbox.ymin()-1 
       << " "
       << scale*bbox.xmax()+1 
       << " " 
       << scale*bbox.ymax()+1
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


void StraightSkeletonFree(int* numFace_i, double* xf, double* yf)
{
  delete [] xf;
  delete [] yf;
  delete [] numFace_i;
}

typedef Sls::Face_iterator   Face_iterator ;
typedef Sls::Face_handle     Face_handle ;
typedef Sls::Halfedge_handle Halfedge_handle ;

int StraightSkeleton( int np
                                          , int* np_i
                                          , double* xp
                                          , double* yp
                                          , int& numFaces
                                          , int& numVertices
                                          , int*& numFace_i
                                          , double*& xf
                                          , double*& yf
                                          , int dumpEPS
                                          )
{
  int result = 0 ;

  numFace_i = NULL ;
  xf = yf = NULL ;

  try
  {
    double scale = 1.0;

    SlsBuilder ssb ;  

    
    int currentPoint = 0;
    for(int i = 0; i < np; i++)
    {
      std::vector<Point> points(np_i[i]);
      for(int j=0; j < np_i[i]; j++)
      {
        Point p(xp[currentPoint], yp[currentPoint]);
        points[j] = p;
        ++currentPoint;
      }

      if( ! CGAL::is_simple_2(points.begin(),points.end()))
      {
        std::cerr << "Polygon " << i << "  is not simple" << std::endl;
        return 0;
      }

      if(  CGAL::orientation_2(points.begin(),points.end()) != ( i == 0 ? CGAL::COUNTERCLOCKWISE 
                                                                        : CGAL::CLOCKWISE 
                                                               ) 
        )
           ssb.enter_contour(points.rbegin(),points.rend());
      else ssb.enter_contour(points.begin(),points.end());
    }  
      

    // Construct the skeleton
    boost::shared_ptr<Sls> ss = ssb.construct_skeleton();
      
    // Proceed only if the skeleton was correctly constructed.
    if ( ss )
    {
        // We first count the points
        numFaces= (int)ss->size_of_faces();
        numFace_i = new int[numFaces];

        numVertices = 0;
        int currentFace = 0;
        for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit)
        {
	         int count = 0;
          Halfedge_handle h = fit->halfedge();
	         Halfedge_handle done;
	         done = h;
	         do
          {
	           count++;
            h = h->next();
	         } while(h != done);
       	  numVertices += count;

	         numFace_i[currentFace] = count;
	         ++currentFace;
        }
          
        // Allocate the x and y array and traverse the faces again
        xf = new double[numVertices];
        yf = new double[numVertices];

        int currentVertex = 0;

        for(Face_iterator fit = ss->faces_begin(); fit != ss->faces_end(); ++fit)
        {
	         Halfedge_handle h = fit->halfedge();
	         Halfedge_handle done;
	         done = h;
	         do
          {
       	    xf[currentVertex] =  h->vertex()->point().x();
	           yf[currentVertex] =  h->vertex()->point().y();
            ++currentVertex;
            h = h->next();
	         }
          while(h != done);
        }

        //
        int vi = 0 ;
        for ( int fi = 0 ; fi < numFaces ; ++ fi )
        {
          double lastx = xf[vi] ;
          double lasty = yf[vi];

          double firstx ;
          double firsty ;

          for ( int fvi = 0 ; fvi < numFace_i[fi] ; ++ fvi )
          {
             firstx = xf[vi];
             firsty = yf[vi];
             ++ vi ;
          }


          Face_iterator ff = ss->faces_begin() ;
          std::advance(ff,fi) ;
          Face_handle fh = ff;
          double _lastx  = fh->halfedge()->vertex()->point().x();
          double _lasty  = fh->halfedge()->vertex()->point().y();
          double _firstx = fh->halfedge()->opposite()->vertex()->point().x();
          double _firsty = fh->halfedge()->opposite()->vertex()->point().y();

          if ( firstx != _firstx || firsty != _firsty || lastx != _lastx || lasty != _lasty )
          {
            std::cout << "face " << fi << " edge: (" << firstx << "," << firsty << ")->(" << lastx << "," << lasty << ") mismatch:" 
                      << " (" << _firstx << "," << _firsty << ")->(" << _lastx << "," << _lasty << ")\n" ;
          }

        }
        //


        result = 1 ;
    }
  }
  catch ( std::exception const& e ) 
  {
    std::cerr << "Exception thrown: " << e.what() << std::endl ;
    StraightSkeletonFree(numFace_i,xf,yf);  
  }
  catch ( ... ) 
  {
    std::cerr << "Unknown exception thrown." << std::endl ;
    StraightSkeletonFree(numFace_i,xf,yf);  
  }

  return result ;
}

void TestDLL ( Region const& aRegion )
{
  int np = 1 ;

  int nv = 0 ;
  for( Region::const_iterator bit = aRegion.begin(), ebit = aRegion.end() ; bit != ebit ; ++ bit )
     nv += (*bit)->size();

  int*    np_i = new int   [np];
  double* xp   = new double[nv];
  double* yp   = new double[nv];

  int  numFaces ;
  int  numVertices ;
  int* numFace_i ;
  double* xf ;
  double* yf ;

  int ci = 0 ;
  int vi = 0 ; 
  for( Region::const_iterator bit = aRegion.begin(), ebit = aRegion.end() ; bit != ebit ; ++ bit )
  {
    np_i[ci] = (*bit)->size();

    for ( Polygon::const_iterator vit = (*bit)->begin() ; vit != (*bit)->end() ; ++ vit )
    {
       xp[vi] = vit->x() ;
       yp[vi] = vit->y() ;
       ++ vi ;
    }
    ++ ci ;
  }

  StraightSkeleton(np,np_i,xp,yp,numFaces,numVertices,numFace_i,xf,yf,0);
  StraightSkeletonFree(numFace_i,xf,yf);

  delete[] np_i ;
  delete[] xp ;
  delete[] yp ;
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
TestDLL(*lRegion);
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
      case 'e' : sDumpEPS = true ; break ;
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
    cout << "USAGE: test_sls_builder <options> folder file0 file1 ... fileN" << endl
         << "  <options>: " << endl
         << "     -e  Dumps result into an .eps file." << endl
         << "     -h       Prints this usage summary." << endl ;
  }

  return sFailed == 0 ? 0 : 1 ;
}
