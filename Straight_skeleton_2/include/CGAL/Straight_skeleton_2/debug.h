// Copyright (c) 2007 Fernando Luis Cacciola Carballal. All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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

#ifndef CGAL_STRAIGHT_SKELETON_DEBUG_H
#define CGAL_STRAIGHT_SKELETON_DEBUG_H 1

#include <CGAL/config.h>

#ifdef CGAL_USE_CORE
#  include <CGAL/CORE_BigFloat.h>
#endif

#if    defined(CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) \
    || defined(CGAL_POLYGON_OFFSET_ENABLE_TRACE) \
    || defined(CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE) \
    || defined(CGAL_STRAIGHT_SKELETON_ENABLE_VALIDITY_TRACE) \
    || defined(CGAL_STRAIGHT_SKELETON_ENABLE_INTRINSIC_TESTING)
#
#  define CGAL_STSKEL_TRACE_ON
#
#  include<string>
#  include<iostream>
#  include<sstream>
#  include<iomanip>
#  define CGAL_STSKEL_TRACE(m) \
     { \
       std::ostringstream ss ; \
       ss << m ; \
       std::string s = ss.str(); \
       Straight_skeleton_external_trace(s); \
     }

#include <boost/optional.hpp>
#include <boost/intrusive_ptr.hpp>

template<class T>
inline std::string o2str( boost::optional<T> const& o )
{
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  if ( o )
       ss << *o ;
  else ss << "·NONE·" ;
  return ss.str();
}

template<class T>
inline std::string ptr2str( boost::intrusive_ptr<T> const& ptr )
{
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  if ( ptr )
       ss << *ptr ;
  else ss << "·NULL·" ;
  return ss.str();
}

template<class N>
inline std::string n2str( N const& n )
{
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  
  ss << CGAL_NTS to_double(n) ;
  
  return ss.str();
}


#if 0 //CGAL_USE_CORE

inline CORE::BigFloat to_big_float( CGAL::MP_Float const& n )
{
  return n.to_rational<CORE::BigFloat>() ;
}

inline CORE::BigFloat to_big_float( CGAL::Quotient<CGAL::MP_Float> const& q )
{
  CORE::BigFloat n = to_big_float(q.numerator  ()) ;
  CORE::BigFloat d = to_big_float(q.denominator()) ;
  if ( !d.isZeroIn())
       return n / d ;
  else return CORE::BigFloat(std::numeric_limits<double>::infinity());
}

template<class NT>
inline CORE::BigFloat to_big_float( NT const& n )
{
  return CORE::BigFloat( CGAL_NTS to_double(n) ) ;
}


inline std::string n2str( CGAL::MP_Float const& n )
{
  std::ostringstream ss ; 
  ss << to_big_float(n) ;
  return ss.str();
}

inline std::string n2str( CGAL::Quotient< CGAL::MP_Float > const& n )
{
  std::ostringstream ss ; 
  ss << to_big_float(n) ;
  return ss.str();
}
#else
inline std::string n2str( CGAL::MP_Float const& n )
{
  std::ostringstream ss ; ss << std::setprecision(19) ;
  ss << CGAL_NTS to_double(n) ;
  return ss.str();
}

inline std::string n2str( CGAL::Quotient< CGAL::MP_Float > const& n )
{
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  ss << CGAL_NTS to_double(n) ;
  return ss.str();
}
#endif

template<class XY>
inline std::string xy2str( XY const& xy )
{
  std::ostringstream ss ; 
  ss << "(" << n2str(xy.x()) << "," << n2str(xy.y()) << ")" ;
  return ss.str();
}
template<class D>
inline std::string dir2str( D const& d )
{
  std::ostringstream ss ; 
  ss << "(" << n2str(d.dx()) << "," << n2str(d.dy()) << ")" ;
  return ss.str();
}
template<class P>
inline std::string p2str( P const& p )
{
  std::ostringstream ss ; 
  ss << "(" << n2str(p.x()) << "," << n2str(p.y()) << ")" ;
  return ss.str();
}
template<class OP>
inline std::string op2str( OP const& op )
{
  return op ? p2str(*op) : std::string("·NONE·");
}
template<class V>
inline std::string v2str( V const& v )
{
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  ss << "V" << v.id() << " " << p2str(v.point()) << " [" << v.time() << "]" ;
  return ss.str();
}
template<class VH>
inline std::string vh2str( VH const& vh )
{
  VH null ;
  return vh != null ? v2str(*vh) : "NULL_VERTEX_HANDLE" ;
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
  std::ostringstream ss ; ss << std::setprecision(19)  ;
  if ( e.is_bisector() )
  {
    ss << "B" << e.id()
       << "[E" << e.defining_contour_edge()->id() 
       << ",E" << e.opposite()->defining_contour_edge()->id() << "]"
       << " (/" << ( e.slope() == CGAL::ZERO ? "·" : ( e.slope() == CGAL::NEGATIVE ? "-" : "+" ) )
       << " " << e.opposite()->vertex()->time() << "->" << e.vertex()->time() << ")" ; 
  }
  else
  {
    ss << "E" << e.id() ;
  }
  ss << " " << s2str(e.opposite()->vertex()->point(),e.vertex()->point()) ;
  return ss.str();
}

template<class EH>
inline std::string eh2str( EH const& eh )
{
  EH null ;
  return eh != null ? e2str(*eh) : "NULL_HALFEDGE_HANDLE" ;
}

template<class BH>
inline std::string newb2str( char const* name, BH const& b )
{
  std::ostringstream ss ; 
  
  ss << "New Bisector " 
     << name 
     << " is B" << b->id()
     << " [E" << b->defining_contour_edge()->id() 
     << ",E" << b->opposite()->defining_contour_edge()->id()
     << "] {B" << b->prev()->id() 
     << "->N"  << b->prev()->vertex()->id() 
     << "->B" << b->id() 
     << "->N" << b->vertex()->id() 
     << "->B" << b->next()->id()
     << "}" ;
     
  return ss.str();
}

template<class VH, class Triedge>
inline std::string newn2str( char const* name, VH const& v, Triedge const& aTriedge )
{
  std::ostringstream ss ; 

  ss << "New Node " << name <<" is N" << v->id() << " at " << v->point()
     << " [E" << aTriedge.e0()->id()
     << ",E" << aTriedge.e1()->id()
     << ",E" << aTriedge.e2()->id()
     << "] incident halfedge: B" << v->halfedge()->id()
     << "  primary bisector: B" << v->primary_bisector()->id() ;
     
  return ss.str();
}

#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
#  define CGAL_STSKEL_DEBUG_CODE(code) code
#  define CGAL_STSKEL_BUILDER_TRACE(l,m) if ( l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#  define CGAL_STSKEL_BUILDER_TRACE_IF(c,l,m) if ( (c) && l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#else
#  define CGAL_STSKEL_DEBUG_CODE(code) 
#  define CGAL_STSKEL_BUILDER_TRACE(l,m)
#  define CGAL_STSKEL_BUILDER_TRACE_IF(c,l,m)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_TRACE
#  define CGAL_POLYOFFSET_DEBUG_CODE(code) code
#  define CGAL_POLYOFFSET_TRACE(l,m) if ( l <= CGAL_POLYGON_OFFSET_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#else
#  define CGAL_POLYOFFSET_DEBUG_CODE(code)
#  define CGAL_POLYOFFSET_TRACE(l,m)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_TRAITS_ENABLE_TRACE
bool sEnableTraitsTrace = false ;
#  define CGAL_STSKEL_TRAITS_ENABLE_TRACE sEnableTraitsTrace = true ;
#  define CGAL_STSKEL_TRAITS_ENABLE_TRACE_IF(cond) if ((cond)) sEnableTraitsTrace = true ;
#  define CGAL_STSKEL_TRAITS_DISABLE_TRACE sEnableTraitsTrace = false;
#  define CGAL_STSKEL_TRAITS_TRACE(m) \
     if ( sEnableTraitsTrace ) \
     { \
       std::ostringstream ss ; \
       ss << m ; \
       std::string s = ss.str(); \
       Straight_skeleton_traits_external_trace(s); \
     }
#else
#  define CGAL_STSKEL_TRAITS_ENABLE_TRACE
#  define CGAL_STSKEL_TRAITS_ENABLE_TRACE_IF(cond)
#  define CGAL_STSKEL_TRAITS_DISABLE_TRACE
#  define CGAL_STSKEL_TRAITS_TRACE(m)
#endif


#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_VALIDITY_TRACE
#  define CGAL_STSKEL_VALIDITY_TRACE(m) CGAL_STSKEL_TRACE(m)
#  define CGAL_STSKEL_VALIDITY_TRACE_IF(cond,m) if ( cond ) CGAL_STSKEL_VALIDITY_TRACE(m)
#else
#  define CGAL_STSKEL_VALIDITY_TRACE(m) 
#  define CGAL_STSKEL_VALIDITY_TRACE_IF(cond,m)
#endif


#ifdef CGAL_STRAIGHT_SKELETON_PROFILING_ENABLED // Reserved use. DO NOT define this macro switch
#  include<string>
#  include<iostream>
#  include<sstream>

namespace CGAL {

namespace CGAL_STRAIGHT_SKELETON_i_profiling
{

template<class NT> char const* kernel_type() { return typeid(NT).name() ; }

template<> char const* kernel_type<double>              () { return "double" ;   }
template<> char const* kernel_type<Interval_nt_advanced>() { return "Interval" ; }
template<> char const* kernel_type< Quotient<MP_Float> >() { return "MP_Float" ; }
template<> char const* kernel_type<CORE::Expr>          () { return "Expr" ;     }

} // CGAL_STRAIGHT_SKELETON_i_profiling

} // end namespace CGAL

#define CGAL_STSKEL_ASSERT_PREDICATE_RESULT(expr,K,pred,error) \
        { \
          std::ostringstream predss ; \
          predss << CGAL_STRAIGHT_SKELETON_i_profiling::kernel_type< typename K::FT >() << " . " << pred ; \
          std::string preds = predss.str(); \
          if ( ! is_certain((expr)) ) \
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

#endif // CGAL_STRAIGHT_SKELETON_DEBUG_H //
// EOF //

 
