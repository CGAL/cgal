// Copyright (c) 2007  GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the soNTware.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Fernando Cacciola <fernando.cacciola@gmail.com>

#include "basics.h"
#include "test_self_intersection.h"


// This is here only to allow a breakpoint to be placed so I can trace back the problem.
void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  cerr << "CGAL error: " << what << " violation!" << endl 
       << "Expr: " << expr << endl
       << "File: " << file << endl 
       << "Line: " << line << endl;
  if ( msg != 0)
    cerr << "Explanation:" << msg << endl;
    
  throw std::runtime_error(expr);  
}

namespace SMS = CGAL::Surface_mesh_simplification ;

typedef SMS::Edge_profile<Surface> Profile ;

class Dump_link_builder : public CGAL::Modifier_base<Surface::HalfedgeDS> 
{
  map<int,int> mMap ;
  
  int mVIdx ;
  
protected:

  typedef boost::graph_traits<Surface> GraphTraits ;
  
  typedef GraphTraits::out_edge_iterator out_edge_iterator ;
  typedef GraphTraits::edge_descriptor   edge_descriptor ;
  typedef GraphTraits::vertex_descriptor vertex_descriptor ;
  
  typedef CGAL::Polyhedron_incremental_builder_3<Surface::HalfedgeDS> Builder ;
  
  Dump_link_builder() : mVIdx(0) {}
  
  void add_vertex( Builder& aBuilder, Vertex_const_handle v )
  {
    aBuilder.add_vertex(v->point());
    mMap.insert( make_pair(v->id(),mVIdx++) );
  }
  
  void add_triangle( Builder& aBuilder, Profile::Triangle const& aTri )
  {
    aBuilder.begin_facet();
    aBuilder.add_vertex_to_facet( mMap[aTri.v0->id()] );
    aBuilder.add_vertex_to_facet( mMap[aTri.v1->id()] );
    aBuilder.add_vertex_to_facet( mMap[aTri.v2->id()] );
    aBuilder.end_facet();
  }
};

class Dump_edge_link_builder : public Dump_link_builder
{ 
  Profile const& mProfile ;
  
  
public:

  Dump_edge_link_builder( Profile const& aProfile ) : mProfile(aProfile) {}
  
  void operator()( Surface::HalfedgeDS& hds) 
  {
    Builder B( hds, true);
    
    B.begin_surface( 2 + mProfile.link().size(), mProfile.triangles().size() );
    
    this->add_vertex(B,mProfile.v0());
    this->add_vertex(B,mProfile.v1());
    
    for ( Profile::vertex_descriptor_vector::const_iterator vit = mProfile.link().begin(); vit != mProfile.link().end(); ++ vit )
      this->add_vertex(B, *vit );
    
    for ( Profile::Triangle_vector::const_iterator tit = mProfile.triangles().begin(); tit != mProfile.triangles().end(); ++ tit )
      this->add_triangle(B,*tit);
    
    B.end_surface();
  }
};

class Dump_vertex_link_builder : public Dump_link_builder
{
  Surface*      mECM ;
  Vertex_handle mV ;
  
  Surface& ecm() { return *mECM ; }
  
public:

  Dump_vertex_link_builder( Surface& aECM, Vertex_handle aV) : mECM(&aECM), mV(aV) {}
  
  void operator()( Surface::HalfedgeDS& hds ) 
  {
    Builder B( hds, true);
    
    B.begin_surface( 1 + boost::out_degree(mV,ecm()), boost::out_degree(mV,ecm()) );
    
    this->add_vertex(B,mV);
    Profile::Triangle_vector triangles ;
    
    out_edge_iterator eb, ee ; 
    for ( tie(eb,ee) = boost::out_edges(mV,ecm()) ; eb != ee ; ++ eb )
    {
      edge_descriptor out_edge1 = *eb ;
      edge_descriptor out_edge2 = out_edge1->opposite()->next(); 
      vertex_descriptor v1 = boost::target(out_edge1,ecm());
      vertex_descriptor v2 = boost::target(out_edge2,ecm());
       
      this->add_vertex(B,v1);
      
      triangles.push_back( Profile::Triangle(mV,v1,v2) ) ;
    }  
       
    for ( Profile::Triangle_vector::const_iterator tit = triangles.begin(); tit != triangles.end(); ++ tit )
      this->add_triangle(B,*tit);
    
    B.end_surface();
  }
};

void dump_edge_link ( Profile const& aProfile, string aName )
{
  Surface lLink ;
  Dump_edge_link_builder lBuilder(aProfile) ;
    
  lLink.delegate(lBuilder);
    
  ofstream out(aName.c_str());
  out << lLink ;    
}

void dump_vertex_link ( Profile const& aProfile, Vertex_handle aV, string aName )
{
  Surface lLink ;
  Dump_vertex_link_builder lBuilder(aProfile.surface(),aV) ;
    
  lLink.delegate(lBuilder);
    
  ofstream out(aName.c_str());
  out << lLink ;    
}

template<class T>
string opt2str ( optional<T> const& o )
{
  ostringstream ss ;
  if ( o )
       ss << *o ; 
  else ss << "<none>" ;
  return ss.str(); 
}

template<class N>
string float2str ( N n )
{
  return boost::str( boost::format("%+08.8e") % n ) ;
}

template<class P>
string point2str ( P const& p )
{
  return boost::str( boost::format( "(%+08.8e %+08.8e %+08.8e)") % p.x() % p.y() % p.z() ) ;
}

template<class P>
string optpoint2str ( optional<P> const& p )
{
  ostringstream ss ;
  if ( p )
       ss << point2str(*p);
  else ss << "<none>" ;
  
  return ss.str(); 
}
template<class N>
string optfloat2str ( optional<N> const& n )
{
  ostringstream ss ;
  if ( n )
       ss << float2str(*n);
  else ss << "<none>" ;
  
  return ss.str(); 
}

template<class V>
string vertex2str ( V const& v )
{
  ostringstream ss ;
  ss << "[V" << v->id() << point2str(v->point()) << "]" ;
  return ss.str(); 
}

template<class E>
string edge2str ( E const& e )
{
  ostringstream ss ;
  ss << "{E" << e->id() << vertex2str(e->opposite()->vertex()) << "->" << vertex2str(e->vertex()) << "}" ;
  return ss.str(); 
}

template<class T> ostream&  operator << ( ostream& os, optional<T> const& o ) { return os << opt2str(o); }

string normalize_EOL ( string line )
{
  string::size_type l = line.length();
  string::size_type d = ( l > 0 && line[l-1] == '\r' ) ? 1 : 0 ; 
  return line.substr(0, l-d ) ;
}

#define REPORT_ERROR(msg) error(__FILE__,__LINE__,0,msg);

#define REPORT_ERROR2(pred,msg) REPORT_ERROR(msg)
          
#define CHECK_MSG(pred,msg) if (!(pred)) REPORT_ERROR2(#pred,msg) 
         
#define CHECK(pred) CHECK_MSG(pred,string(""))
        
#define CHECK_EQUAL(x,y)       CHECK_MSG(((x)==(y)), str(format("Assertion failed: %1%(=%2%)==%3%(=%4%)") % (#x) % (x) % (#y) % (y) ) )
#define CHECK_NOT_EQUAL(x,y)   CHECK_MSG(((x)!=(y)), str(format("Assertion failed: %1%(=%2%)!=%3%(=%4%)") % (#x) % (x) % (#y) % (y) ) )

class Visitor
{
  
public :
  
  void OnStarted( Surface& ) 
  { 
    mStep = 0 ; 
  } 
  
  void OnFinished ( Surface& aSurface )
  { 
    CHECK(aSurface.is_valid());
  } 
  
  void OnStopConditionReached( Profile const& ) {} 
  
  void OnCollected( Profile const& aProfile, optional<NT> const& aCost )
  {
    TRACE( str ( format("Collecting %1% : cost=%2%") % edge2str(aProfile.v0_v1()) % optfloat2str(aCost) ) ) ;
  }                
  
  void OnSelected( Profile const&, optional<NT> const&, size_t, size_t ) 
  {
    ::CGALi::cgal_enable_ecms_trace = mStep == 0 ;
  }
  
  void OnCollapsing( Profile const& aProfile, optional<Point> const& aP ) 
  {
    TRACE( str ( format("S %1% - Collapsing %2% : placement=%3%") % mStep % edge2str(aProfile.v0_v1()) % optpoint2str(aP) ) ) ;
//    if ( mStep == 11 )
//      dump_edge_link(aProfile,"step-11-before.off");
  }
  
  void OnCollapsed( Profile const& aProfile, Vertex_handle aV )
  {
    if ( Is_self_intersecting( aProfile.surface() ) )
    {
      dump_vertex_link( aProfile, aV, str( format( "step-%1%-after.off" ) % mStep ) ) ;
      cerr << "Resulting surface self-intersects." << endl ;
    }  
    
    ++ mStep ;
  }

  void OnNonCollapsable( Profile const& ) {}                

private:

  void error ( char const* file, int line, char const* pred, string msg )
  {
    cerr << "ERROR in " << file << " at " << line << endl ;
    if ( pred )
      cerr << "  Assertion failed: " << pred << endl ;
    cerr << "  " << msg << endl ;
   
     
    throw runtime_error("");
  }
  
private:
    
  int mStep ;   
} ;

bool sWriteResult = false ;
bool sUseMP       = false ; 
int  sStop        = 1 ;

bool Test ( string aName )
{
  bool rSucceeded = false ;
  
  try
  {
    string off_name = aName ;
    
    ifstream off_is(off_name.c_str());
    if ( off_is )
    {
      Surface lSurface; 
      off_is >> lSurface ;
      if ( lSurface.is_valid() )
      {
        if ( lSurface.is_pure_triangle() )
        {
          cerr << "Processing " << aName << " (" << ( lSurface.size_of_halfedges() / 2 ) << " edges)" << endl ;
          
          Visitor vis ;
          
          set_halfedgeds_items_id(lSurface);
          
          SMS::Count_stop_predicate<Surface> stop(sStop);

          Real_timer t ; t.start();    

          if ( sUseMP )
          {
            SMS::Edge_length_cost  <Surface> cost ;
            SMS::Midpoint_placement<Surface> placement ;
            
            edge_collapse(lSurface,stop,get_cost(cost).get_placement(placement).visitor(&vis) );
          }
          else
          {
            SMS::LindstromTurk_cost     <Surface> cost ;
            SMS::LindstromTurk_placement<Surface> placement ;
            
            edge_collapse(lSurface,stop,get_cost(cost).get_placement(placement).visitor(&vis) );
          } 

          t.stop();
                       
          rSucceeded = lSurface.is_valid() && lSurface.is_pure_triangle() ;
          
          cerr << "\r" << aName << ( rSucceeded ? " succeeded" : "FAILED" ) << ". Elapsed time=" << t.time() << " seconds." << endl ;

          if ( sWriteResult )
          {
            string out_name = off_name ;
            out_name.replace( off_name.find_last_of('.'), 15, ".simplified.off");
            ofstream out(out_name.c_str());
            out << lSurface ;
            cerr << "Result written into: " << out_name << endl ;
          }
        }
        else
        {
          cerr << "Surface is not triangulated (has faces with more than 3 sides): " << aName << endl ;
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
  }
  catch ( std::exception const& x ) 
  {
    string what(x.what());
    if ( what.length() > 0 )
      cerr << "Exception caught: " << what << endl ;
  }
  
  
  return rSucceeded ;
}

int main( int argc, char** argv ) 
{
  set_error_handler  (error_handler);
  set_warning_handler(error_handler);
  
  cout << setprecision(4);
  cerr << setprecision(4);
  
  vector<string> lCases ;
  
  for ( int i = 1 ; i < argc ; ++i )
  {
    
    if ( argv[i][0] == '-' )
    {
      switch(argv[i][1])
      { 
        case 'w': sWriteResult = true ; break ;
        case 's': sStop = atoi(&argv[i][2]) ; break ;
        case 'm': sUseMP = true ; break ;
      }
    }
    else
    {
      string c( normalize_EOL( string(argv[i]) ) ) ;
      string::size_type pos = c.find_last_of(".") ;
      if ( pos != string::npos )
      {
        string ext = c.substr(pos);
        if ( ext == ".off" )
          lCases.push_back(c);
      }
    }
  }
   
  if ( lCases.size() == 0 )
  {
    cout << "collapse_edge_test file0 [-w] [-sNUM_EDGES] [-m] file1 ... fileN\n" ;
    return 1 ;
  } 
  else
  {
    unsigned lOK = 0 ;
    for ( vector<string>::const_iterator it = lCases.begin(); it != lCases.end() ; ++ it )
    {
     if ( Test(*it) )
       ++ lOK ;
    }  
      
    cout << endl
         << lOK                    << " cases succedded." << endl
         << (lCases.size() - lOK ) << " cases failed." << endl ;
            
    return lOK == lCases.size() ? 0 : 1 ;
  }
}

#ifdef WHY_DO_THESE_DONT_EXIST_IN_CGAL_UNDER_LINUX
namespace CGAL
{

void assertion_fail( const char* expr, const char* file, int line )
{
  assertion_fail(expr,file,line,"");
}
void precondition_fail( const char* expr, const char* file, int line )
{
  precondition_fail(expr,file,line,"");
}
void postcondition_fail( const char* expr, const char* file, int line )
{
  postcondition_fail(expr,file,line,"");
}

}

#endif

// EOF //


