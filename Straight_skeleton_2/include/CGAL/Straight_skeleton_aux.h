// Copyright (c) 2005, 2006 Fernando Luis Cacciola Carballal. All rights reserved.
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

#if defined(CGAL_STRAIGHT_SKELETON_ENABLE_TRACE) || defined(CGAL_POLYGON_OFFSET_ENABLE_TRACE)
#define CGAL_SLS_ENABLE_TRACE
#endif

#if   defined(CGAL_STRAIGHT_SKELETON_ENABLE_SHOW) \
   || defined(CGAL_POLYGON_OFFSET_ENABLE_SHOW) \
   || defined(CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX)  \
   || defined(CGAL_POLYGON_OFFSET_ENABLE_SHOW_AUX)
#define CGAL_SLS_ENABLE_SHOW
#endif


#ifdef CGAL_SLS_ENABLE_TRACE
#  include<string>
#  include<iostream>
#  include<sstream>
#  define CGAL_SLS_TRACE(m) \
     { \
       std::ostringstream ss ; ss << m ; std::string s = ss.str(); \
       Straight_skeleton_external_trace(s); \
     }
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_TRACE
#  define CGAL_SSBUILDER_TRACE(m) CGAL_SLS_TRACE(m)
#else
#  define CGAL_SSBUILDER_TRACE(m)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW
#  define CGAL_SSBUILDER_SHOW(code) code
#else
#  define CGAL_SSBUILDER_SHOW(code)
#endif

#ifdef CGAL_STRAIGHT_SKELETON_ENABLE_SHOW_AUX
#  define CGAL_SSBUILDER_SHOW_AUX(code) code
#else
#  define CGAL_SSBUILDER_SHOW_AUX(code)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_TRACE
#  define CGAL_POLYOFFSET_TRACE(m) CGAL_SLS_TRACE(m)
#else
#  define CGAL_POLYOFFSET_TRACE(m)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_SHOW
#  define CGAL_POLYOFFSET_SHOW(code) code
#else
#  define CGAL_POLYOFFSET_SHOW(code)
#endif

#ifdef CGAL_POLYGON_OFFSET_ENABLE_SHOW_AUX
#  define CGAL_POLYOFFSET_SHOW_AUX(code) code
#else
#  define CGAL_POLYOFFSET_SHOW_AUX(code)
#endif


CGAL_BEGIN_NAMESPACE

#ifdef CGAL_SLS_ENABLE_SHOW
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

}
#endif

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

#ifdef CGAL_SLS_PROFILING_ENABLED // Reserved use. DO NOT define this macro switch
#  include<string>
#  include<iostream>
#  include<sstream>

namespace CGAL_SLS_i_profiling
{

template<class NT> char const* kernel_type() { return typeid(NT).name() ; }

template<> char const* kernel_type<double>              () { return "double" ;   }
template<> char const* kernel_type<Interval_nt_advanced>() { return "Interval" ; }
template<> char const* kernel_type< Quotient<MP_Float> >() { return "MP_Float" ; }
template<> char const* kernel_type<CORE::Expr>          () { return "Expr" ;     }

}

//
// Undefined identifier: register_predicate_failure ??
// Then you mistakenly defined the macro switch CGAL_SLS_PROFILING_ENABLED, turn it off.
//
#define CGAL_SLS_ASSERT_PREDICATE_RESULT(expr,K,pred,error) \
        { \
          std::ostringstream predss ; \
          predss << CGAL_SLS_i_profiling::kernel_type< typename K::FT >() << " . " << pred ; \
          std::string preds = predss.str(); \
          if ( is_indeterminate((expr)) ) \
          { \
            std::ostringstream errss  ; errss << error ; std::string errs = errss.str(); \
            register_predicate_failure(preds,errs); \
          } \
          else register_predicate_success(preds); \
        }
#else
#define CGAL_SLS_ASSERT_PREDICATE_RESULT(expr,K,pred,error)
#endif

CGAL_END_NAMESPACE

namespace boost
{
inline void intrusive_ptr_add_ref( CGAL::Ref_counted_base const* p ) { p->AddRef(); }
inline void intrusive_ptr_release( CGAL::Ref_counted_base const* p ) { p->Release(); }
} // namespace boost


#undef CGAL_SLS_ENABLE_TRACE
#undef CGAL_SLS_ENABLE_SHOW

#endif // CGAL_STRAIGHT_SKELETON_AUX_H //
// EOF //

