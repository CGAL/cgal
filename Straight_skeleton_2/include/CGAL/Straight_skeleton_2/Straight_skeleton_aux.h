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

enum Triedge_collinearity
{ 
    TRIEDGE_COLLINEARITY_NONE
  , TRIEDGE_COLLINEARITY_01
  , TRIEDGE_COLLINEARITY_12
  , TRIEDGE_COLLINEARITY_02
  , TRIEDGE_COLLINEARITY_ALL
} ;

namespace CGALi 
{

template <>
struct Minmax_traits< Triedge_collinearity >
{
  static const Triedge_collinearity min = TRIEDGE_COLLINEARITY_NONE;
  static const Triedge_collinearity max = TRIEDGE_COLLINEARITY_ALL;
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
#define CGAL_STSKEL_ENABLE_TRACE
#endif

#if   defined(CGAL_STRAIGHT_SKELETON_ENABLE_SHOW) \
   || defined(CGAL_POLYGON_OFFSET_ENABLE_SHOW) \
   || defined(CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX)  \
   || defined(CGAL_POLYGON_OFFSET_ENABLE_SHOW_AUX)
#define CGAL_STSKEL_ENABLE_SHOW
#endif


#ifdef CGAL_STSKEL_ENABLE_TRACE
#  include<string>
#  include<iostream>
#  include<sstream>
#  define CGAL_STSKEL_TRACE(m) \
     { \
       std::ostringstream ss ; ss << m ; std::string s = ss.str(); \
       Straight_skeleton_external_trace(s); \
     }

#  define CGAL_STSKEL_DEBUG_CODE(code) code
#else
#  define CGAL_STSKEL_DEBUG_CODE(code) 
#endif

#ifdef CGAL_STSKEL_ENABLE_TRACE
#  define CGAL_STSKEL_BUILDER_TRACE(l,m) if ( l <= CGAL_STRAIGHT_SKELETON_ENABLE_TRACE ) CGAL_STSKEL_TRACE(m)
#else
#  define CGAL_STSKEL_BUILDER_TRACE(l,m)
#endif

#ifdef CGAL_STSKEL_ENABLE_SHOW
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

#ifdef CGAL_STSKEL_ENABLE_SHOW

CGAL_BEGIN_NAMESPACE

namespace SS_IO_AUX
{
  class ScopedDrawing
  {
    public :

      virtual ~ScopedDrawing()
      {
        if ( mID != -1 )
          Straight_skeleton_external_undraw_object(mID) ;
      }

      void Release() { mID = -1 ; }
    protected :

      ScopedDrawing ( int aID ) : mID(aID) {}

    private :

      int mID ;
  } ;

  class ScopedPointDrawing : public ScopedDrawing
  {
    public :

    template<class Point_2>
    ScopedPointDrawing( Point_2 const& aP, CGAL::Color aColor, char const* aLayer )
      :
      ScopedDrawing
      (
        Straight_skeleton_external_draw_point(  to_double( aP.x() )
                                               ,to_double( aP.y() )
                                               ,aColor
                                               ,aLayer
                                             )
      )
    {}
  } ;

  class ScopedSegmentDrawing : public ScopedDrawing
  {
    public :

    template<class Point_2>
    ScopedSegmentDrawing( Point_2 const& aS, Point_2 const& aT, CGAL::Color aColor, char const* aLayer )
      :
      ScopedDrawing
      (
        Straight_skeleton_external_draw_segment(  to_double( aS.x() )
                                                 ,to_double( aS.y() )
                                                 ,to_double( aT.x() )
                                                 ,to_double( aT.y() )
                                                 ,aColor
                                                 ,aLayer
                                               )
      )
    {}
  } ;

} // namespace SS_IO_AUX

CGAL_END_NAMESPACE

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
#undef CGAL_STSKEL_ENABLE_SHOW

#endif // CGAL_STRAIGHT_SKELETON_AUX_H //
// EOF //

