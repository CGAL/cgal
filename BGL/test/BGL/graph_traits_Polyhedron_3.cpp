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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Straight_skeleton_2/test/Straight_skeleton_2/test_sls_builder.cpp $
// $Id: test_sls_builder.cpp 32700 2006-07-24 22:31:02Z fcacciola $
//
// Author(s)     : Fernando Cacciola <fernando_cacciola@ciudad.com.ar>
//

#include <cstdio> 
#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Simple_cartesian<double> Kernel ;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron ;

int sOK     = 0 ; 
int sFailed = 0 ;

using namespace std ;
 
struct Match
{
  std::size_t num_vertices ;
};

istream& operator >> ( istream& is, Match& rMatch )
{
  return is ;
}

#define CHECK_EQUAL(x,y) \
        if ( (x) != (y) ) \
        { \
          cerr << "Assertion failure: " << (x) << "==" << (y) << endl \
               << "File:" << __FILE__ << endl \
               << "Line:" << __LINE__ << endl ; \
          throw 0 ; \
        }

void test ( Polyhedron const& aPoly, Match const& aMatch )
{
  using namespace boost ;
  
  CHECK_EQUAL( num_vertices(aPoly), aPoly.size_of_vertices () ) 
  CHECK_EQUAL( num_edges   (aPoly), aPoly.size_of_halfedges() ) 
}

bool test( string off_file )
{
  bool rContinue = true ;
  
  std::size_t extpos = off_file.find_last_of('.') ;
  if ( extpos != string::npos && off_file.substr(extpos) == ".off" )
  {
    ifstream is(off_file.c_str());
    if ( is )
    {
      string traits_file = off_file.substr(0,extpos)+".traits";
      ifstream match(traits_file.c_str());
      if ( match )
      {
        Polyhedron lPoly ;
        is >> lPoly ;
        
        Match lMatch ;
        match >> lMatch ;
        
        bool ok = false ;
        try
        {
          test(lPoly,lMatch) ;
          ok = true ;
        } catch(...) {}  
        
        if ( ok )
             ++ sOK ;
        else ++ sFailed ;  
        
        cout << ( ok ? "OK" : "FAILED!" ) << endl ;
      }
      else cerr << "Unable to load .traits file: " << traits_file << endl ;
    }
    else cerr << "Unable to load input .off file: " << off_file << endl ;
  }
  else cerr << "Input file must have .off extension: " << off_file << endl ;
  
  return rContinue ;
}

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!\n"
            << "Expr: " << expr << endl
            << "File: " << file << endl 
            << "Line: " << line << endl;
  if ( msg != 0)
      cerr << "Explanation:" << msg << endl;
      
  // Avoid an abort()    
  throw runtime_error(msg);     
}

int main( int argc, char const* argv[] )
{
  cout << "Testing Polyhedron's Graph Traits\n";
  
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  bool print_usage = false ;
  bool nop = false ;
  
  string folder = "" ;
  
  for ( int i = 1 ; i < argc ; ++ i )  
  {
    if ( argv[i][0] == '#' )
    {
      folder = string(&argv[i][1]);
      cout << "Input folder: " << folder << endl ;
      break ;
    }  
  }
  
  vector<string> samples ;
  
  for ( int i = 1 ; i < argc ; ++ i )  
    if ( argv[i][0] != '-' && argv[i][0] != '#' )
      samples.push_back(folder+ string(argv[i])); 
  
  if ( samples.size() > 0 ) 
  {
    for ( vector<string>::const_iterator it = samples.begin() ; it != samples.end() ; ++ it )
    {
      if ( !nop )
      {
        if (!test(*it) )
          break ;
      }
      else
        cout << *it << endl ; 
    }    

    int lTotal = sOK + sFailed ;
     
    if ( lTotal > 0 )
    {
      cout << "Total cases: " << lTotal << endl
                << "Succeeded cases: " << sOK << endl
                << "Failed cases: " << sFailed << endl
                << "Failure ratio: " << ((double)sFailed/lTotal*100.0) << "%\n" ;
    } 
    
  }
  else print_usage = true ;

  if ( print_usage )
  {
    cout << "USAGE: graph_traits_Polyhedron_3 #folder file0 file1 ... fileN\n"
              << "  If file is '*' then all files with extension .off in the current folder are loaded\n" ;
  }

  return sFailed == 0 ? 0 : 1 ;
}

