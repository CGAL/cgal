#include <CGAL/assertions_behaviour.h>

bool test( string off_file )
{
  bool rContinue = true ;
  
  std::size_t extpos = off_file.find_last_of('.') ;
  if ( extpos != string::npos && off_file.substr(extpos) == ".off" )
  {
    ifstream is(off_file.c_str());
    if ( is )
    {
      Polyhedron lPoly ;
      is >> lPoly ;
      
      set_halfedgeds_items_id(lPoly);
      
      bool ok = test(lPoly) ;
      
      if ( ok )
           ++ sOK ;
      else ++ sFailed ;  
      
      cout << ( ok ? "OK" : "FAILED!" ) << endl ;
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

int aux_main( int argc, char const* argv[] )
{
  cout << "Testing " << TEST_NAME << endl ;
  
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
    cout << "USAGE: " << TEST_PROGRAM << " graph_traits_Polyhedron_3 #folder file0 file1 ... fileN\n"
              << "  If file is '*' then all files with extension .off in the current folder are loaded\n" ;
  }

  return sFailed == 0 ? 0 : 1 ;
}

