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
// $URL$
// $Id$
// 
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>


#include <iostream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/tokenizer.hpp>

#define CGAL_CHECK_EXPENSIVE

#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_BGL.h>
#include <CGAL/Polyhedron_extended_BGL.h>
#include <CGAL/Polyhedron_BGL_properties.h>
#include <CGAL/IO/Polyhedron_geomview_ostream.h>
#include <CGAL/IO/Polyhedron_iostream.h>

//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_LT_TRACE 4
//#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_TRACE 3
#define CGAL_SURFACE_SIMPLIFICATION_ENABLE_AUDIT

void Surface_simplification_external_trace( std::string s )
{
  static std::ofstream lout("tsms_log.txt");
  lout << s << std::flush << std::endl ;
  std::printf("%s\n",s.c_str());
} 

int exit_code = 0 ;

#include <CGAL/Surface_mesh_simplification_vertex_pair_collapse.h>

#include <CGAL/Surface_mesh_simplification/Policies/Construct_minimal_collapse_data.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_length_cost.h>
#include <CGAL/Surface_mesh_simplification/Policies/Midpoint_vertex_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/LindstromTurk.h>
#include <CGAL/Surface_mesh_simplification/Policies/Count_stop_pred.h>

using namespace std ;
using namespace boost ;
using namespace CGAL ;

//typedef Simple_homogeneous<int>                        Kernel;
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

CGAL_END_NAMESPACE  

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


void ParseAudit ( string name )
{
  ifstream is(name.c_str());
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
  else cerr << "Warning: Audit file " << name << " doesn't exist." << endl ;
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

void CGAL_TSMS_audit( Vertex_handle    const& p 
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

// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!" << endl 
       << "Expr: " << expr << endl
       << "File: " << file << endl 
       << "Line: " << line << endl;
  if ( msg != 0)
    cerr << "Explanation:" << msg << endl;
}

using namespace CGAL::Triangulated_surface_mesh::Simplification ;

char const* matched_alpha ( bool matched )
{
  return matched ? "matched" : "UNMATCHED" ; 
}

bool Test ( string name )
{
  bool rSucceeded = false ;
  
  string off_name    = name+string(".off");
  string audit_name  = name+string(".audit");
  string result_name = name+string(".out.off");
  
  ifstream off_is(off_name.c_str());
  if ( off_is )
  {
    Polyhedron lP; 
    
    off_is >> lP ;
    
    cout << "Testing Lindstrom Turk simplification of surace with " << (lP.size_of_halfedges()/2) << " edges..." << endl ;
  
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

    sAuditData.clear();
    //ParseAudit(audit_name);
      
    cout << "Audit data loaded." << endl ;
    
    typedef LindstromTurk_collapse_data<Polyhedron> Collapse_data ;
    
    Construct_LindstromTurk_collapse_data<Collapse_data> Construct_collapse_data ;
    LindstromTurk_cost                   <Collapse_data> Get_cost ;
    LindstromTurk_vertex_placement       <Collapse_data> Get_vertex_point ;
    Count_stop_condition                 <Collapse_data> Should_stop(0);
        
    Collapse_data::Params lParams;
    
    sAuditReport.clear();
    int r = vertex_pair_collapse(lP,Construct_collapse_data,&lParams,Get_cost,Get_vertex_point,Should_stop);
            
    ofstream off_out(result_name.c_str(),ios::trunc);
    off_out << lP ;
    
    cout << "Finished...\n"
         << r << " edges removed.\n"
         << lP.size_of_vertices() << " vertices.\n"
         << (lP.size_of_halfedges()/2) << " edges.\n"
         << lP.size_of_facets() << " triangles.\n" 
         << ( lP.is_valid() ? " valid\n" : " INVALID!!\n" ) ;

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
  }
  else
  {
    cerr << "Unable to open test file " << name << endl ;
  }              
  
  return rSucceeded ;
}

int main( int argc, char** argv ) 
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);
  
  vector<string> cases ;
    
  for ( int i = 1 ; i < argc ; ++i )
    cases.push_back( string(argv[i]) ) ;
    
  if ( cases.size() == 0 )
    cases.push_back( string("data/sample5") ) ;
  
  unsigned lOK = 0 ;
  for ( vector<string>::const_iterator it = cases.begin(); it != cases.end() ; ++ it )
  {
    if ( Test(*it) )
      ++ lOK ;
  }  
  
  cout << endl
       << lOK                   << " cases succedded." << endl
       << (cases.size() - lOK ) << " cases failed." << endl ;
        
  return lOK == cases.size() ? 0 : 1 ;
}

// EOF //
