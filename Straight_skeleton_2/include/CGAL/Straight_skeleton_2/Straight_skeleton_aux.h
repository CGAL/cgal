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
#ifndef CGAL_STRAIGHT_SKELETON_AUX_H
#define CGAL_STRAIGHT_SKELETON_AUX_H 1

#include <boost/optional/optional.hpp>
#include <boost/none.hpp>

//
// The heap objects used in this implementation are intrusively reference counted. Thus, they inherit from Ref_counted_base.
//
CGAL_BEGIN_NAMESPACE

namespace CGAL_SS_i
{

//
// This record encapsulates the defining contour halfedges for a node (both contour and skeleton)
//
template<class Handle_>
struct Triedge
{
  typedef Handle_ Handle ;
  
  typedef Triedge<Handle> Self ;
  
  Triedge() {}

  // Contour nodes (input polygon vertices) have just 2 defining contour edges    
  Triedge ( Handle aE0, Handle aE1 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
  }              
  
  // Skeleton nodes (offset polygon vertices) have 3 defining contour edges    
  Triedge ( Handle aE0, Handle aE1 , Handle aE2 )
  {
    mE[0] = aE0 ;
    mE[1] = aE1 ;
    mE[2] = aE2 ;
  }              
  
  Handle e( unsigned idx ) const { CGAL_assertion(idx<3); return mE[idx]; }
  
  Handle e0() const { return e(0); }
  Handle e1() const { return e(1); }
  Handle e2() const { return e(2); }
  
  bool is_valid() const { return e0() != e1() && e1() != e2() ; }
  
  // returns 1 if aE is one of the halfedges stored in this triedge, 0 otherwise.
  int contains ( Handle aE ) const
  {
    return aE == e0() || aE == e1() || aE == e2() ? 1 : 0 ;
  }

  // Returns the number of common halfedges in the two triedges x and y
  static int CountInCommon( Self const& x, Self const& y )
  {
    return x.contains(y.e0()) + x.contains(y.e1()) + x.contains(y.e2()) ; 
  }
  
  // Returns true if the triedges store the same 3 halfedges (in any order)
  friend bool operator == ( Self const& x, Self const& y ) { return CountInCommon(x,y) == 3 ; }
  
  friend bool operator != ( Self const& x, Self const& y ) { return !(x==y) ; }
  
  friend Self operator & ( Self const& x, Self const& y )
  {
    return Self(x.e0(), x.e1(), ( x.e0() == y.e0() || x.e1() == y.e0() ) ? y.e1() : y.e0() ) ;
  }
  
  friend std::ostream& operator<< ( std::ostream& ss, Self const& t )
  {
    return ss << "{E" << t.e0()->id() << ",E" << t.e1()->id() << ",E" << t.e2()->id() << "}" ;
  }
  
  Handle mE[3];
} ;

} // namespace CGAL_SS_i

enum Trisegment_collinearity
{ 
    TRISEGMENT_COLLINEARITY_NONE
  , TRISEGMENT_COLLINEARITY_01
  , TRISEGMENT_COLLINEARITY_12
  , TRISEGMENT_COLLINEARITY_02
  , TRISEGMENT_COLLINEARITY_ALL
} ;

static char const* trisegment_collinearity_to_string( Trisegment_collinearity c )
{
  switch ( c )
  {
    case TRISEGMENT_COLLINEARITY_NONE : return "<>" ;  
    case TRISEGMENT_COLLINEARITY_01   : return "<0,1>" ; 
    case TRISEGMENT_COLLINEARITY_12   : return "<1,2>" ; 
    case TRISEGMENT_COLLINEARITY_02   : return "<0,2>" ; 
    case TRISEGMENT_COLLINEARITY_ALL  : return "<0,1,2>" ; 
  }
  
  return "!!UNKNOWN COLLINEARITY!!" ;
}
namespace CGALi 
{

template <>
struct Minmax_traits< Trisegment_collinearity >
{
  static const Trisegment_collinearity min = TRISEGMENT_COLLINEARITY_NONE;
  static const Trisegment_collinearity max = TRISEGMENT_COLLINEARITY_ALL;
};

}

class Ref_counted_base
{
private:
  mutable long mCount ;
  Ref_counted_base( Ref_counted_base const &);
  Ref_counted_base& operator=( Ref_counted_base const &);
protected:
  Ref_counted_base(): mCount(0) {}
  virtual ~Ref_counted_base() {}
public:
    void AddRef() const { ++mCount; }
    void Release() const
      {
        if( --mCount == 0 )
          delete this;
      }
};

CGAL_END_NAMESPACE

namespace boost
{
inline void intrusive_ptr_add_ref( CGAL::Ref_counted_base const* p ) { p->AddRef(); }
inline void intrusive_ptr_release( CGAL::Ref_counted_base const* p ) { p->Release(); }
} // namespace boost

//
// The rest of this header contains tracing, debugging and profiling stuff.

#if defined(CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) || defined(CGAL_POLYGON_OFFSET_ENABLE_TRACE)
#  include<string>
#  include<iostream>
#  include<sstream>
#  define CGAL_STSKEL_TRACE(m) \
     { \
       std::ostringstream ss ; \
       ss << m ; \
       std::string s = ss.str(); \
       Straight_skeleton_external_trace(s); \
     }

template<class N>
inline std::string n2str( N const& n )
{
  std::ostringstream ss ; 
  ss << CGAL_NTS to_double(n) ;
  return ss.str();
}
template<class P>
inline std::string p2str( P const& p )
{
  std::ostringstream ss ; 
  ss << "(" << n2str(p.x()) << "," << n2str(p.y()) << ")" ;
  return ss.str();
}
template<class V>
inline std::string v2str( V const& v )
{
  std::ostringstream ss ; 
  ss << "V" << v.id() << " " << p2str(v.point());
  return ss.str();
}
template<class P>
inline std::string s2str( P const& s, P const& t )
{
  std::ostringstream ss ; 
  ss << "{" << p2str(s) << "-" << p2str(t) << "}" ;
  return ss.str();
}
template<class S>
inline std::string s2str( S const& seg ) { return s2str(seg.source(),seg.target()); }

template<class E>
inline std::string e2str( E const& e )
{
  std::ostringstream ss ; 
  if ( e.is_bisector() )
  {
    ss << "B" << e.id()
       << "[E" << e.defining_contour_edge()->id() 
       << ",E" << e.opposite()->defining_contour_edge()->id() << "]";
  }
  else
  {
    ss << "E" << e.id() ;
  }
  ss << " " << s2str(e.vertex()->point(),e.opposite()->vertex()->point()) ;
  return ss.str();
}
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
#  define CGAL_STSKEL_DEBUG_CODE(code) code
#  define CGAL_STSKEL_BUILDER_TRACE(l,m) if ( l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#else
#  define CGAL_STSKEL_BUILDER_TRACE(l,m)
#  define CGAL_STSKEL_DEBUG_CODE(code) 
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
#  define CGAL_STSKEL_BUILDER_SHOW(code) { code }
#else
#  define CGAL_STSKEL_BUILDER_SHOW(code)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_TRACE
#  define CGAL_POLYOFFSET_TRACE(l,m) if ( l <= CGAL_POLYGON_OFFSET_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#else
#  define CGAL_POLYOFFSET_TRACE(l,m)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_SHOW
#  define CGAL_POLYOFFSET_SHOW(code) { code }
#else
#  define CGAL_POLYOFFSET_SHOW(code)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_PROFILING_ENABLED // Reserved use. DO NOT define this macro switch
#  include<string>
#  include<iostream>
#  include<sstream>

CGAL_BEGIN_NAMESPACE

namespace CGAL_STRAIGHT_SKELETON_i_profiling
{

template<class NT> char const* kernel_type() { return typeid(NT).name() ; }

template<> char const* kernel_type<double>              () { return "double" ;   }
template<> char const* kernel_type<Interval_nt_advanced>() { return "Interval" ; }
template<> char const* kernel_type< Quotient<MP_Float> >() { return "MP_Float" ; }
template<> char const* kernel_type<CORE::Expr>          () { return "Expr" ;     }

} // CGAL_STRAIGHT_SKELETON_i_profiling

CGAL_END_NAMESPACE

#define CGAL_STSKEL_ASSERT_PREDICATE_RESULT(expr,K,pred,error) \
        { \
          std::ostringstream predss ; \
          predss << CGAL_STRAIGHT_SKELETON_i_profiling::kernel_type< typename K::FT >() << " . " << pred ; \
          std::string preds = predss.str(); \
          if ( is_indeterminate((expr)) ) \
          { \
            std::ostringstream errss  ; errss << error ; std::string errs = errss.str(); \
            register_predicate_failure(preds,errs); \
          } \
          else register_predicate_success(preds); \
        }

#define CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(expr,K,cons,error) \
        { \
          std::ostringstream consss ; \
          consss << CGAL_STRAIGHT_SKELETON_i_profiling::kernel_type< typename K::FT >() << " . " << cons ; \
          std::string conss = consss.str(); \
          if ( !(expr) ) \
          { \
            std::ostringstream errss  ; errss << error ; std::string errs = errss.str(); \
            register_construction_failure(conss,errs); \
          } \
          else register_construction_success(conss); \
        }
#else

#define CGAL_STSKEL_ASSERT_PREDICATE_RESULT(expr,K,pred,error)
#define CGAL_STSKEL_ASSERT_CONSTRUCTION_RESULT(expr,K,cons,error)

#endif

#undef CGAL_STSKEL_ENABLE_TRACE

#endif // CGAL_STRAIGHT_SKELETON_AUX_H //
// EOF //

