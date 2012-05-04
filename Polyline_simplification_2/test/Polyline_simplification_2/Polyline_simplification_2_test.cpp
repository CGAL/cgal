#pragma warning( disable : 4503 )

#include <fstream>
#include <vector>

// CGAL headers
#include <CGAL/Bbox_2.h>
#include <CGAL/assertions_behaviour.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>

#define CGAL_TESTING_POLYLINE_SIMPLIFICATION
#define CGAL_POLYLINE_SIMPLIFICATION_TRACE_LEVEL 15

bool lAppToLog = false ;
void Polyline_simplification_2_external_trace( std::string m )
{
  std::ofstream out("polysim_log.txt", ( lAppToLog ? std::ios::app | std::ios::ate : std::ios::trunc | std::ios::ate ) );
  out << std::setprecision(19) << m << std::endl << std::flush ; 
  lAppToLog = true ;
}

#define STR(m) static_cast<std::ostringstream&>( std::ostringstream() << std::string() << m << std::endl ).str()

int sErrors = 0 ;

void error ( std::string msg )
{
  std::cerr << msg << std::endl ;
  ++ sErrors ;
}


void error_handler ( char const* what, char const* expr, char const* file, int line, char const* msg )
{
  error( STR( "CGAL error: " << what << " violation!" << std::endl
            << "Expr: " << expr << std::endl
            << "File: " << file << std::endl
            << "Line: " << line << std::endl
            << "Explanation:" << ( msg ? msg : "" )
            )
       ) ;       
       
  throw std::runtime_error("CGAL ERROR");     
}


#include <CGAL/Polyline_simplification_2/Stop_below_count_threshold.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>
#include <CGAL/Polyline_simplification_2.h>

using namespace std ;

namespace PS2 = CGAL::Polyline_simplification_2 ;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef K::Point_2   Point_2;

typedef CGAL::Simplify_polylines_2<K> PS;

typedef PS::Vertex_handle          Vertex_handle ;
typedef PS::Open_polyline          Open_polyline ;
typedef PS::Open_polyline_handle   Open_polyline_handle ;
typedef PS::Closed_polyline        Closed_polyline ;
typedef PS::Closed_polyline_handle Closed_polyline_handle ;

void dump_polygon( Closed_polyline_handle aPoly )
{
  Closed_polyline::Circulator h = aPoly->circulator();
  Closed_polyline::Circulator c = h ;
  
  do
  {
    std::cout << "V" << c->id() << " " ;
  }  
  while ( ++ c != h ) ;
  
  std::cout << std::endl ;
}

typedef std::vector<int>     ID_vector ;
typedef std::vector<Point_2> Point_vector ;

Point_vector dump_vertex_point_sequence( Open_polyline_handle aPoly )
{
  Point_vector rSeq ;

  for ( Open_polyline::Iterator it = aPoly->begin() ; it != aPoly->end() ; ++ it )
  {
    rSeq.push_back( it->point() ) ;
  } 
  
  return rSeq ;
}

Point_vector dump_vertex_point_sequence( Closed_polyline_handle aPoly )
{
  Point_vector rSeq ;

  Closed_polyline::Circulator h = aPoly->circulator();
  Closed_polyline::Circulator c = h ;
  
  do
  {
    rSeq.push_back( c->point() ) ;
  }  
  while ( ++ c != h ) ;
  
  
  return rSeq ;
}

ID_vector dump_vertex_id_sequence( Open_polyline_handle aPoly )
{
  ID_vector rSeq ;

  for ( Open_polyline::Iterator it = aPoly->begin() ; it != aPoly->end() ; ++ it )
  {
    rSeq.push_back( it->id() ) ;
  } 
  
  return rSeq ;
}

ID_vector dump_vertex_id_sequence( Closed_polyline_handle aPoly )
{
  ID_vector rSeq ;

  Closed_polyline::Circulator h = aPoly->circulator();
  Closed_polyline::Circulator c = h ;
  
  do
  {
    rSeq.push_back( c->id() ) ;
  }  
  while ( ++ c != h ) ;
  
  
  return rSeq ;
}

template<class T>
std::string to_string( T const& t )
{
  std::ostringstream ss ;
  ss << t ;
  return ss.str();  
}

std::string to_string( Vertex_handle const& v )
{
  std::ostringstream ss ;
  ss << "[" << v->id() << "@" << v->point() << "]" ;
  return ss.str();  
}

template<class It>
std::string sequence_to_string( It beg, It end )
{
  std::ostringstream ss ;
  ss << "[" ;
  while ( beg != end )
    ss << to_string(*(beg++)) << "," ;
  ss << "]";
  return ss.str();  
}

template<class InputIteratorA, class InputIteratorB>
void check_sequence ( InputIteratorA got_beg, InputIteratorA got_end, InputIteratorB expected_beg, InputIteratorB expected_end )
{
  std::size_t got_s      = std::distance(got_beg, got_end) ;
  std::size_t expected_s = std::distance(expected_beg, expected_end) ;
  
  bool ok = false ;
  
  if ( got_s == expected_s )
  {
    std::pair<InputIteratorA,InputIteratorB> mr = std::mismatch(got_beg,got_end,expected_beg);
    
    if ( mr.first == got_end || mr.second == expected_end )
      ok = true ;
  }
  
  if ( !ok )
    error( STR("Vertex sequences mismatch:\n  Got=" << sequence_to_string(got_beg, got_end) << "\n  Expected=" << sequence_to_string(expected_beg, expected_end)) );   
} 

void test_crosses_case_0()
{
  std::cout << "Testing crossing polygons case 0" << std::endl ;
  
  PS ps;

  Point_2 lPtsA[4] = { Point_2(0,-10)
                     , Point_2(10,-10)  
                     , Point_2(10,20)  
                     , Point_2(0,20)  
                     } ;
     
  Point_2 lPtsB[4] = { Point_2(-10,0)
                     , Point_2(20,0)
                     , Point_2(20,10)
                     , Point_2(-10,10)
                     } ;
                  
  Closed_polyline_handle lPolyA = ps.insert_polygon(lPtsA,lPtsA+4) ;
  Closed_polyline_handle lPolyB = ps.insert_polygon(lPtsB,lPtsB+4) ;

  ps.simplify( PS2::Stop_below_count_threshold(0), PS2::Squared_distance_cost()) ;

  ID_vector lSeqA = dump_vertex_id_sequence(lPolyA);
  ID_vector lSeqB = dump_vertex_id_sequence(lPolyB);
  
  int lExpectedA[] = {9,10,11,8} ;
  int lExpectedB[] = {8,9,10,11};
  
  check_sequence(lSeqA.begin(), lSeqA.end(), lExpectedA, lExpectedA+4);
  check_sequence(lSeqB.begin(), lSeqB.end(), lExpectedB, lExpectedB+4);
  
  Point_2 lResA[4] = { Point_2(10,0)
                     , Point_2(10,10)  
                     , Point_2(0,10)  
                     , Point_2(0,0)  
                     } ;
     
  Point_2 lResB[4] = { Point_2(0,0)
                     , Point_2(10,0)
                     , Point_2(10,10)
                     , Point_2(0,10)
                     } ;
  Point_vector lPA = dump_vertex_point_sequence(lPolyA);
  Point_vector lPB = dump_vertex_point_sequence(lPolyB);
  
  check_sequence(lPA.begin(), lPA.end(), lResA, lResA+4);
  check_sequence(lPB.begin(), lPB.end(), lResB, lResB+4);
    
}

void test_crosses_case_1()
{
  std::cout << "Testing crossing polygons case 1" << std::endl ;
  
  PS ps;

  Point_2 lPtsA[4] = { Point_2(0,0)
                     , Point_2(10,0)  
                     , Point_2(10,10)  
                     , Point_2(0,10)  
                     } ;
     
  Point_2 lPtsB[4] = { Point_2(5,-5)
                     , Point_2(15,-5)
                     , Point_2(15,15)
                     , Point_2(5,15)
                     } ;
                  
  Closed_polyline_handle lPolyA = ps.insert_polygon(lPtsA,lPtsA+4) ;
  Closed_polyline_handle lPolyB = ps.insert_polygon(lPtsB,lPtsB+4) ;

  ps.simplify( PS2::Stop_below_count_threshold(0), PS2::Squared_distance_cost()) ;

  ID_vector lSeqA = dump_vertex_id_sequence(lPolyA);
  ID_vector lSeqB = dump_vertex_id_sequence(lPolyB);
  
  int lExpectedA[] = {9,8,3} ;
  int lExpectedB[] = {5,6,8,9};
  
  check_sequence(lSeqA.begin(), lSeqA.end(), lExpectedA, lExpectedA+3);
  check_sequence(lSeqB.begin(), lSeqB.end(), lExpectedB, lExpectedB+4);
  
  Point_2 lResA[3] = { Point_2(5,0)
                     , Point_2(5,10)  
                     , Point_2(0,10)  
                     } ;
     
  Point_2 lResB[4] = { Point_2(15,-5)
                     , Point_2(15,15)
                     , Point_2(5,10)
                     , Point_2(5,0)
                     } ;
                     
  Point_vector lPA = dump_vertex_point_sequence(lPolyA);
  Point_vector lPB = dump_vertex_point_sequence(lPolyB);
  
  check_sequence(lPA.begin(), lPA.end(), lResA, lResA+3);
  check_sequence(lPB.begin(), lPB.end(), lResB, lResB+4);
  
}

void test_crosses_case_2()
{
  std::cout << "Testing crossing polylines case 0" << std::endl ;
  
  PS ps;

  Point_2 lPtsA[4] = { Point_2(0,0)
                     , Point_2(10,10)  
                     , Point_2(20,0)  
                     } ;
  Point_2 lPtsB[4] = { Point_2(10,20)
                     , Point_2(10,0)  
                     } ;
                  
  Open_polyline_handle lPolyA = ps.insert_polyline(lPtsA,lPtsA+3) ;
  Open_polyline_handle lPolyB = ps.insert_polyline(lPtsB,lPtsB+2) ;

  ps.simplify( PS2::Stop_below_count_threshold(0), PS2::Squared_distance_cost()) ;

  ID_vector lSeqA = dump_vertex_id_sequence(lPolyA);
  ID_vector lSeqB = dump_vertex_id_sequence(lPolyB);
  
  int lExpectedA[] = {0,1,2} ;
  int lExpectedB[] = {3,1,4} ;
  
  check_sequence(lSeqA.begin(), lSeqA.end(), lExpectedA, lExpectedA+3);
  check_sequence(lSeqB.begin(), lSeqB.end(), lExpectedB, lExpectedB+3);

  Point_2 lResA[3] = { Point_2(0,0)
                     , Point_2(10,10)  
                     , Point_2(20,0)  
                     } ;
     
  Point_2 lResB[3] = { Point_2(10,20)
                     , Point_2(10,10)
                     , Point_2(10,0)
                     } ;
                     
  Point_vector lPA = dump_vertex_point_sequence(lPolyA);
  Point_vector lPB = dump_vertex_point_sequence(lPolyB);
  
  check_sequence(lPA.begin(), lPA.end(), lResA, lResA+3);
  check_sequence(lPB.begin(), lPB.end(), lResB, lResB+3);
  
}

struct Visitor
{
  Visitor( ID_vector* aSeq ) : mSeq(aSeq) {}
  
  void OnStarted() const {} 
  
  void OnFinished() const {} 
  
  void OnStopConditionReached() const {} 
  
  template<class VH>
  void OnCollected( VH const& ) const {}                
  
  template<class VH>
  void OnSelected( VH const&, boost::optional<double> const&, unsigned, unsigned ) const {}                
  
  template<class VH>
  void OnRemoving( VH const&, VH const& aV, VH const& ) const { mSeq->push_back(aV->id()) ; }          
  
  template<class VH>
  void OnRemoved( VH const&, VH const&) const {}        
  
  template<class VH>
  void OnNonRemovable( VH const& ) const {}                
  
  ID_vector* mSeq ;
  
} ;

void test_removal()
{
  std::cout << "Testing removal sequence" << std::endl ;
  
  PS ps;

  Point_2 lPts[10] = { Point_2(0,0)
                     , Point_2(1,1)  
                     , Point_2(3,-3)  
                     , Point_2(9,9)  
                     , Point_2(0,18)  
                     } ;
                  
  Open_polyline_handle lPoly = ps.insert_polyline(lPts,lPts+5) ;
    
  ID_vector lRemovedSeq ;
  
  Visitor vis(&lRemovedSeq) ;
  
  ps.simplify( PS2::Stop_below_count_threshold(7), PS2::Squared_distance_cost(), vis) ;
  
  int lExpected[] = {1,2,3} ;
  
  check_sequence(lRemovedSeq.begin(), lRemovedSeq.end(), lExpected, lExpected+3);
  
  Point_2 lRes[7] = { Point_2(0,0)
                    , Point_2(0,18)  
                    } ;

  Point_vector lPP = dump_vertex_point_sequence(lPoly);
  
  check_sequence(lPP.begin(), lPP.end(), lRes, lRes+2);

}

void test_fixed()
{
  std::cout << "Testing fixed vertices" << std::endl ;
  
  PS ps;

  Point_2 lPts[10] = { Point_2(0,0)
                     , Point_2(1,1)  
                     , Point_2(2,0)  
                     , Point_2(6,0)  
                     , Point_2(7,2)  
                     , Point_2(8,0)  
                     , Point_2(13,0)  
                     , Point_2(14,4)  
                     , Point_2(15,0)  
                     , Point_2(19,0)  
                     } ;
                  
  Open_polyline_handle lPoly = ps.insert_polyline(lPts,lPts+10) ;
    
  for ( Open_polyline::Iterator it = lPoly->begin() ; it != lPoly->end() ; ++ it )
    if ( it->point().y() == 0 )
      it->set_is_fixed(true);
  
  ID_vector lRemovedSeq ;
  
  Visitor vis(&lRemovedSeq) ;
  
  ps.simplify( PS2::Stop_below_count_threshold(7), PS2::Squared_distance_cost(), vis) ;
  
  int lExpected[] = {1,4,7} ;
  
  check_sequence(lRemovedSeq.begin(), lRemovedSeq.end(), lExpected, lExpected+3);
  
  Point_2 lRes[7] = { Point_2(0,0)
                    , Point_2(2,0)  
                    , Point_2(6,0)  
                    , Point_2(8,0)  
                    , Point_2(13,0)  
                    , Point_2(15,0)  
                    , Point_2(19,0)  
                    } ;

  Point_vector lPP = dump_vertex_point_sequence(lPoly);
  
  check_sequence(lPP.begin(), lPP.end(), lRes, lRes+7);

}

int main(int argc, char **argv)
{
  int rOK = 1 ;
  
  CGAL::set_error_handler  (error_handler);
  CGAL::set_warning_handler(error_handler);
  
  try 
  { 
    test_crosses_case_0();
    test_crosses_case_1();
    test_crosses_case_2();
    test_removal();
    test_fixed();

    rOK = sErrors == 0  ;
  }
  catch(exception& e) 
  {
    cerr << e.what() << "\n";
  }
  
  return rOK ;
}
