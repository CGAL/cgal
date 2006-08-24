// // Copyright (c) 2002  Max Planck Institut fuer Informatik (Germany).
// All rights reserved.
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
// $URL: svn+ssh://fcacciola@scm.gforge.inria.fr/svn/cgal/trunk/Surface_mesh_simplification/test/Surface_mesh_simplification/LT_edge_collapse_test.cpp $
// $Id: LT_edge_collapse_test.cpp 32177 2006-07-03 11:55:13Z fcacciola $
// 
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>


#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Real_timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_BGL.h>
#include <CGAL/Polyhedron_extended_BGL.h>
#include <CGAL/Polyhedron_BGL_properties.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Constrained_triangulation_2.h>

//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE 4
//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 4

#define TRACK_STATS
#define SHOW_STATS
//#define AUDIT

void Surface_simplification_external_trace( std::string s )
{
  static std::ofstream lout("tsms_log.txt");
  lout << s << std::flush << std::endl ;
  std::printf("%s\n",s.c_str());
} 

int exit_code = 0 ;

#include <CGAL/Surface_mesh_simplification_edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Polyhedron_is_vertex_fixed_map.h>
#include <CGAL/Surface_mesh_simplification/Polyhedron_edge_cached_pointer_map.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/LindstromTurk.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_pred.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_pred.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

typedef Simple_cartesian<double> Kernel;
typedef Kernel::Vector_3         Vector;
typedef Kernel::Point_3          Point;

template <class Refs, class Traits>
struct My_vertex : public HalfedgeDS_vertex_base<Refs,Tag_true,Point> 
{
  typedef HalfedgeDS_vertex_base<Refs,Tag_true,Point> Base ;
 
  My_vertex() : ID(-1), IsFixed(false) {} 
  My_vertex( Point p ) : Base(p), ID(-1), IsFixed(false) {}

  int id() const { return ID ; }
    
  int  ID; 
  bool IsFixed ;
} ;

template <class Refs, class Traits>
struct My_halfedge : public HalfedgeDS_halfedge_base<Refs> 
{ 
  My_halfedge() 
   : 
     ID(-1) 
  {}
 
  int id() const { return ID ; }
  
  int ID; 
  
  void* cached_pointer ;
};

template <class Refs, class Traits>
struct My_face : public HalfedgeDS_face_base<Refs,Tag_true,typename Traits::Plane_3>
{
  typedef HalfedgeDS_face_base<Refs,Tag_true,typename Traits::Plane_3> Base ;
  
  My_face() : ID(-1) {}
  My_face( typename Traits::Plane_3 plane ) : Base(plane) {}
  
  int id() const { return ID ; }
  
  int ID; 
};

struct My_items : public Polyhedron_items_3 
{
    template < class Refs, class Traits>
    struct Vertex_wrapper {
        typedef My_vertex<Refs,Traits> Vertex;
    };
    template < class Refs, class Traits>
    struct Halfedge_wrapper { 
        typedef My_halfedge<Refs,Traits>  Halfedge;
    };
    template < class Refs, class Traits>
    struct Face_wrapper {
        typedef My_face<Refs,Traits> Face;
    };
};

typedef Polyhedron_3<Kernel,My_items> Polyhedron; 

typedef Polyhedron::Vertex                                   Vertex;
typedef Polyhedron::Vertex_iterator                          Vertex_iterator;
typedef Polyhedron::Vertex_handle                            Vertex_handle;
typedef Polyhedron::Vertex_const_handle                      Vertex_const_handle;
typedef Polyhedron::Halfedge_handle                          Halfedge_handle;
typedef Polyhedron::Halfedge_const_handle                    Halfedge_const_handle;
typedef Polyhedron::Edge_iterator                            Edge_iterator;
typedef Polyhedron::Facet_iterator                           Facet_iterator;
typedef Polyhedron::Halfedge_around_vertex_const_circulator  HV_circulator;
typedef Polyhedron::Halfedge_around_facet_circulator         HF_circulator;

CGAL_BEGIN_NAMESPACE

template<>
struct External_polyhedron_get_is_vertex_fixed<Polyhedron>
{
  bool operator() ( Polyhedron const&, Vertex_const_handle v ) const { return v->IsFixed; }
}  ;

template<>
struct External_polyhedron_access_edge_cached_pointer<Polyhedron>
{
  void*& operator() ( Polyhedron&, Halfedge_handle he ) const { return he->cached_pointer; }
}  ;

CGAL_END_NAMESPACE  


#ifdef AUDIT

typedef CGAL::Constrained_triangulation_2<> CT ;
CT ct ;

double sCostMatchThreshold   = 1e-2 ;
double sVertexMatchThreshold = 1e-2 ;

struct Audit_data
{
  Audit_data( Point p, Point q ) : P(p), Q(q) {}
  
  Point            P, Q ;
  optional<double> Cost ;
  optional<Point>  NewVertexPoint ;
} ;
typedef shared_ptr<Audit_data> Audit_data_ptr ;
typedef vector<Audit_data_ptr> Audit_data_vector ;

Audit_data_vector sAuditData ;

Audit_data_ptr lookup_audit ( Point p, Point q )
{ 
  
  for ( Audit_data_vector::iterator it = sAuditData.begin() ; it != sAuditData.end() ; ++ it )
  {
    Audit_data_ptr data = *it ;
    if ( data->P == p && data->Q == q )
      return data ;
  }
  return Audit_data_ptr() ;  
}

Audit_data_ptr get_or_create_audit ( Point p, Point q )
{ 
  Audit_data_ptr data = lookup_audit(p,q);
  if ( !data )
  {
    data = Audit_data_ptr(new Audit_data(p,q));
    sAuditData.push_back(data);
    ct.insert_constriant(p,q));
  }  
  return data ;
}

struct Audit_report
{
  Audit_report( int eid ) 
   : 
     HalfedgeID(eid)
   , CostMatched(false)
   , OrderMatched(false)
   , NewVertexMatched(false)
  {}
 
  
  int              HalfedgeID ;
  Audit_data_ptr   AuditData ; 
  optional<double> Cost ;
  optional<Point>  NewVertexPoint ;
  bool             CostMatched ;
  bool             OrderMatched ;
  bool             NewVertexMatched ; 
} ;

typedef shared_ptr<Audit_report> Audit_report_ptr ;

typedef map<int,Audit_report_ptr> Audit_report_map ;

Audit_report_map sAuditReport ;

typedef char_separator<char> Separator ;
typedef tokenizer<Separator> Tokenizer ;

Point ParsePoint( vector<string> tokens, int aIdx )
{
  double x = atof(tokens[aIdx+0].c_str());
  double y = atof(tokens[aIdx+1].c_str());
  double z = atof(tokens[aIdx+2].c_str());
  
  return Point(x,y,z);
}

void ParseAuditLine( string s )
{
  Tokenizer tk(s,Separator(","));
  vector<string> tokens(tk.begin(),tk.end());
  if ( tokens.size() >= 7 )
  {
    Point p = ParsePoint(tokens,0);
    Point q = ParsePoint(tokens,3);
    
    Audit_data_ptr lData = get_or_create_audit(p,q);
    
    double cost = atof(tokens[6].c_str()) ;
    
    if ( cost != std::numeric_limits<double>::max() )
      lData->Cost = cost ;
    
    if ( tokens.size() == 10 )
      lData->NewVertexPoint = ParsePoint(tokens,7);
  }
  else
    cerr << "WARNING: Invalid audit line: [" << s << "]" << endl ;
}


void ParseAudit ( string aName )
{
  ifstream is(aName.c_str());
  if ( is )
  {
    while ( is )
    {
      string line ;
      getline(is,line);
      if ( line.length() > 0 )
        ParseAuditLine(line);
    }
  }
  else cerr << "Warning: Audit file " << aName << " doesn't exist." << endl ;
}

string to_string( optional<double> const& c )
{
  if ( c )
       return str( format("%1%") % (*c) ) ;
  else return "NONE" ;  
}

string to_string( optional<Point> const& p )
{
  if ( p )
       return str( format("(%1%,%2%,%3%)") % p->x() % p->y() % p->z() ) ;
  else return "NONE" ;  
}

void register_collected_edge( Vertex_handle    const& p 
                            , Vertex_handle    const& q
                            , Halfedge_handle  const& e
                            , optional<double> const& cost
                            , optional<Point>  const& newv
                            )
{
  Audit_report_ptr lReport( new Audit_report(e->ID) ) ;
  
  lReport->Cost           = cost ;
  lReport->NewVertexPoint = newv ;
  
  Audit_data_ptr lData ; //= lookup_audit(p->point(),q->point());
  if ( lData )
  {
    lReport->AuditData = lData ;
    
    if ( !!lData->Cost && !!cost )
    {
      double d = fabs(*lData->Cost - *cost) ;
      if ( d < sCostMatchThreshold )
        lReport->CostMatched = true ;
    }
    else if (!lData->Cost && !cost )
      lReport->CostMatched = true ;
      
    if ( !!lData->NewVertexPoint && !!newv )
    {
      double d = std::sqrt(squared_distance(*lData->NewVertexPoint,*newv));
      
      if ( d < sVertexMatchThreshold )
        lReport->NewVertexMatched = true ;
    }
    else if ( !lData->NewVertexPoint && !newv )
      lReport->NewVertexMatched = true ;
  }
  
  sAuditReport.insert( make_pair(e->ID,lReport) ) ;
} 
#endif

#ifdef TRACK_STATS
int sInitial ;
int sCollected ;
int sProcessed ;
int sCollapsed ;
int sNonCollapsable ;
int sCostUncomputable ;
int sFixed ;
int sRemoved ;
#endif


struct Visitor
{
  void OnStarted( Polyhedron& aSurface ) 
  {
#ifdef TRACK_STATS
    sInitial = aSurface.size_of_halfedges() / 2 ;
#endif
  } 
  
  void OnFinished ( Polyhedron& aSurface ) 
  {
#ifdef TRACK_STATS
    cerr << "\n";
#endif
  } 
  
  void OnStopConditionReached( Polyhedron& aSurface ) {} 
  
  void OnCollected( Halfedge_handle const& aEdge
                  , bool                   aIsFixed
                  , Polyhedron&            aSurface
                  )
  {
#ifdef AUDIT  
    register_collected_edge(aEdge,aCost,aNewVertexPoint); 
#endif

#ifdef TRACK_STATS
    ++sCollected ;
    cerr << "\rEdges collected: " << sCollected;
    if ( aIsFixed )
      ++ sFixed ;
#endif
 }                
  
  void OnProcessed(Halfedge_handle const& aEdge
                  ,Polyhedron&            aSurface
                  ,optional<double>       aCost
                  ,Vertex_handle const&   aVertex
                  )
  {
#ifdef TRACK_STATS
    if ( sProcessed == 0 )
      cerr << "\n";  
    ++ sProcessed ;
    if ( aVertex == Vertex_handle() )
    {
      if ( !aCost )
           ++ sCostUncomputable ;
      else ++ sNonCollapsable ;
    }
    else 
    {
      ++ sCollapsed;
      sRemoved += 3 ;
      cerr << "\r" << ((int)(100.0*((double)sRemoved/(double)sInitial))) << "%" ;
    }
#endif  
  }                
  
} ;

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!" << endl 
       << "Expr: " << expr << endl
       << "File: " << file << endl 
       << "Line: " << line << endl;
  if ( msg != 0)
    cerr << "Explanation:" << msg << endl;
    
  throw std::logic_error(what);  
}

using namespace CGAL::Triangulated_surface_mesh::Simplification::Edge_collapse ;

char const* matched_alpha ( bool matched )
{
  return matched ? "matched" : "UNMATCHED" ; 
}

enum Method { Midpoint, LT_cached, LT_uncached } ;

char const* method_to_string( Method aMethod )
{
  switch(aMethod)
  {
    case Midpoint:  return "midpoint" ; break ;
    case LT_cached: return "LindstromTurk (cached)" ; break ;
    case LT_uncached: return "LindstromTurk (uncached)" ; break ;
  }
  
  return "<unknown>" ;
}

typedef Set_empty_collapse_data             <Polyhedron> P_set_empty_collapse_data ;
typedef Set_full_collapse_data_LindstromTurk<Polyhedron> P_set_full_collapse_data_LT ;

typedef Empty_collapse_data<Polyhedron> P_empty_collapse_data ;
typedef Full_collapse_data<Polyhedron>  P_full_collapse_data ;

typedef LindstromTurk_params LT_params ;
typedef char                 Dummy_params ;

typedef Edge_length_cost<P_empty_collapse_data>   MP_cost ;
typedef LindstromTurk_cost<P_empty_collapse_data> LT_uncached_cost ; 
typedef LindstromTurk_cost<P_full_collapse_data>  LT_cached_cost ;
          
typedef Midpoint_placement     <P_empty_collapse_data> MP_placement ;
typedef LindstromTurk_placement<P_empty_collapse_data> LT_uncached_placement ;
typedef LindstromTurk_placement<P_full_collapse_data>  LT_cached_placement ;

bool Test ( int aStopA, int aStopR, bool aJustPrintSurfaceData, string aName, Method aMethod )
{
  bool rSucceeded = false ;
  
  string off_name    = aName ;
  string audit_name  = aName+string(".audit");
  string result_name = aName+string(".out.off");
  
  ifstream off_is(off_name.c_str());
  if ( off_is )
  {
    Polyhedron lP; 
    
    scan_OFF(off_is,lP,true);
    
    if ( lP.is_valid() )
    {
      if ( lP.is_pure_triangle() )
      {
        if ( !aJustPrintSurfaceData )
        {
          int lFinalEdgesCount ;
          if ( aStopA != -1 )
               lFinalEdgesCount = aStopA ;
          else lFinalEdgesCount = lP.size_of_halfedges() * aStopR / 200 ;
          
          cout << "Testing simplification of surface " << off_name << " using " << method_to_string(aMethod) << "method."  << endl ;
             
          cout << lP.size_of_facets() << " triangles." << endl 
               << (lP.size_of_halfedges()/2) << " edges." << endl 
               << lP.size_of_vertices() << " vertices." << endl 
               << (lP.is_closed() ? "Closed." : "Open." ) << endl 
               << "Final edge count: " << lFinalEdgesCount << endl ;
               
          cout << setprecision(19) ;
        
          int lVertexID = 0 ;
          for ( Polyhedron::Vertex_iterator vi = lP.vertices_begin(); vi != lP.vertices_end() ; ++ vi )
            vi->ID = lVertexID ++ ;    
          
          int lHalfedgeID = 0 ;
          for ( Polyhedron::Halfedge_iterator hi = lP.halfedges_begin(); hi != lP.halfedges_end() ; ++ hi )
            hi->ID = lHalfedgeID++ ;    
          
          int lFacetID = 0 ;
          for ( Polyhedron::Facet_iterator fi = lP.facets_begin(); fi != lP.facets_end() ; ++ fi )
            fi->ID = lFacetID ++ ;    
      
#ifdef AUDIT
          sAuditData .clear();
          sAuditReport.clear();
          ParseAudit(audit_name);
          cout << "Audit data loaded." << endl ;
#endif
          P_set_empty_collapse_data   set_empty_collapse_data ;
          P_set_full_collapse_data_LT set_full_collapse_data_LT ;
          
          MP_cost          get_MP_cost;
          LT_uncached_cost get_LT_uncached_cost;
          LT_cached_cost   get_LT_cached_cost;
          
          MP_placement          get_MP_placement;
          LT_uncached_placement get_LT_uncached_placement;
          LT_cached_placement   get_LT_cached_placement;

                    
          Count_stop_condition<Polyhedron> should_stop(lFinalEdgesCount);
              
          Dummy_params lDummy_params;
          LT_params    lLT_params ; 
          
          Visitor lVisitor ;
      
          int r = -1 ;
          
          Real_timer t ; t.start();    
          switch( aMethod )
          {
            case Midpoint:  
              r = edge_collapse(lP,&lDummy_params,set_empty_collapse_data,get_MP_cost,get_MP_placement,should_stop,&lVisitor);
              break ;
            case LT_cached: 
              r = edge_collapse(lP,&lLT_params,set_full_collapse_data_LT,get_LT_cached_cost,get_LT_cached_placement,should_stop,&lVisitor);
              break ;
            case LT_uncached:
              r = edge_collapse(lP,&lLT_params,set_empty_collapse_data,get_LT_uncached_cost,get_LT_uncached_placement,should_stop,&lVisitor);
              break ;
          }
          t.stop();
                  
          ofstream off_out(result_name.c_str(),ios::trunc);
          off_out << lP ;
          
          cout << "\nFinished...\n"
               << "Ellapsed time: " << t.time() << " seconds.\n" 
               << r << " edges removed.\n"
               << endl
               << lP.size_of_vertices() << " final vertices.\n"
               << (lP.size_of_halfedges()/2) << " final edges.\n"
               << lP.size_of_facets() << " final triangles.\n" 
               << ( lP.is_valid() ? " valid\n" : " INVALID!!\n" ) ;
      
#ifdef SHOW_STATS
          cout << "\n"
               << sProcessed        << " edges processed.\n"
               << sCollapsed        << " edges collapsed.\n" 
               << sNonCollapsable   << " non-collapsable edges.\n"
               << sCostUncomputable << " non-computable edges.\n"
               << sFixed            << " fixed edges.\n"
               << sRemoved          << " edges removed.\n" 
               << (sRemoved/3)      << " vertices removed." ;
#endif
      
#ifdef AUDIT
          unsigned lMatches = 0 ;
                        
          cout << "Audit report:\n" ;
          for ( Audit_report_map::const_iterator ri = sAuditReport.begin() ; ri != sAuditReport.end() ; ++ ri )
          {
            Audit_report_ptr lReport = ri->second ;
            
            if ( lReport->AuditData )
            {
              if ( lReport->NewVertexMatched )
                ++ lMatches ;
                 
              cout << "Collapsed Halfedge " << lReport->HalfedgeID << endl 
                   << "  Cost: Actual=" << to_string(lReport->Cost) << ", Expected=" << to_string(lReport->AuditData->Cost)
                                        << ". " << matched_alpha(lReport->CostMatched) << endl
                   << "  New vertex point: Actual=" << to_string(lReport->NewVertexPoint) << ", Expected=" 
                                                    << to_string(lReport->AuditData->NewVertexPoint) 
                                                    << ". " << matched_alpha(lReport->NewVertexMatched)
                   << endl ;
            }
            else
            {
              cout << "No audit data for Halfedge " << lReport->HalfedgeID << endl ;
                ++ lMatches ;
            }  
          }
          
          rSucceeded = ( lMatches == sAuditReport.size() ) ;
#else
          rSucceeded = true ;
#endif
        }   
        else
        {
          cout << off_name << ": " << lP.size_of_facets() << " triangles, " 
               << (lP.size_of_halfedges()/2) << " edges, "
               << lP.size_of_vertices() << " vertices, "
               << (lP.is_closed() ? "Closed." : "Open." ) << endl ;
               
          rSucceeded = true ;
        }  
      }
      else
      {
        cerr << "Surfaces is not triangulated (has faces with more than 3 sides): " << aName << endl ;
      }
    }
    else
    {
      cerr << "Invalid surface: " << aName << endl ;
    }
  }
  else
  {
    cerr << "Unable to open test file " << aName << endl ;
  }              
  
  return rSucceeded ;
}

bool sPrintUsage = false ;

void add_case( string aCase, vector<string>& rCases )
{
  if ( aCase.find(".") == string::npos )
    aCase += ".off" ;
  
  if ( aCase.find(".off") != string::npos )
  {
    rCases.push_back(aCase);
  }
  else 
  {
    sPrintUsage = true ;  
    cerr << "Invalid input file. Only .off files are supported: " << aCase << endl ;
  }    
}

int main( int argc, char** argv ) 
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);
  
  bool   lJustPrintSurfaceData = false ;
  int    lStopA = -1 ;
  int    lStopR = 20 ;
  Method lMethod = Midpoint ;
  string lFolder =""; 
  vector<string> lCases ;
        
  for ( int i = 1 ; i < argc ; ++i )
  {
    string opt(argv[i]);
    if ( opt[0] == '-' )
    {
      switch(opt[1])
      {
        case 'd' : lFolder = opt.substr(2); break ;
        case 'a' : lStopA = lexical_cast<int>(opt.substr(2)); break;
        case 'r' : lStopR = lexical_cast<int>(opt.substr(2)); break;
        case 'n' : lJustPrintSurfaceData = true ; break ;
        case 'm' : lMethod = (Method)lexical_cast<int>(opt.substr(2)); break ;
        
        default: 
          cerr << "Invalid option: " << opt << endl ;
          sPrintUsage = true ; 
      }
    }
  }
  
  for ( int i = 1 ; i < argc ; ++i )
  {
    string opt(argv[i]);
    if ( opt[0] == '@' )
    {
      string rspname = opt.substr(1) ;
      ifstream rsp(rspname.c_str());
      if ( rsp )
      {
        string line ;
        while ( getline(rsp,line) )
          add_case(lFolder+line,lCases);
      }
      else
      {
        cerr << "Cannot open response file: " << rspname << endl ;
        sPrintUsage = true ; 
      }
    }
    else if ( opt[0] != '-' )
    {
      add_case(lFolder+opt,lCases);
    }
  }
   
  if ( lCases.size() == 0 )
    sPrintUsage = true ;
    
  if ( sPrintUsage )
  {
    cout << "collapse_edge_test <options> file0 file1 ... fileN" << endl 
         << "  options: " << endl
         << "    -m method                    method: 0=midpoint[default] 1=LindstromTurk (cached) 2=LindstromTurk (uncached)" << endl 
         << "    -d folder                    Specifies the folder where the files are located. " << endl 
         << "    -a absolute_max_edge_count   Sets the final number of edges as absolute number." << endl
         << "    -r relative_max_edge_count   Sets the final number of edges as a percentage." << endl
         << "    -n                           Do not simplify but simply report data of surfaces." << endl ;
    
    return 1 ;
  } 
  else
  {
    unsigned lOK = 0 ;
    for ( vector<string>::const_iterator it = lCases.begin(); it != lCases.end() ; ++ it )
    {
     if ( Test( lStopA, lStopR, lJustPrintSurfaceData, *it, lMethod) )
       ++ lOK ;
    }  
      
    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK ) << " cases failed." << endl ;
            
    return lOK == lCases.size() ? 0 : 1 ;
  }
}

// EOF //
