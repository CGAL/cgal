#line 260 "pierce_traits.awi"
#line 18 "code_formatting.awi"
// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision $
// release_date  : $CGAL_Date $
//
// file          : Pierce_rectangles_2_traits.h
// chapter       : $CGAL_Chapter: Geometric Optimisation $
// package       : $CGAL_Package: Matrix_search $
// source        : pcenter.aw
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Hoffmann <hoffmann@inf.ethz.ch>
//
// coordinator   : ETH Zurich (Bernd Gaertner <gaertner@inf.ethz.ch>)
//
// 2-4-Piercing Axis-Parallel 2D-Rectangles
// ============================================================================

#line 264 "pierce_traits.awi"
#line 54 "code_formatting.awi"
#if ! (PIERCE_RECTANGLES_2_TRAITS_H)
#define PIERCE_RECTANGLES_2_TRAITS_H 1

#line 31 "pierce_traits.awi"
#ifndef CGAL_POINT_2_H
#include <CGAL/Point_2.h>
#endif // CGAL_POINT_2_H
#ifndef CGAL_ISO_RECTANGLE_2_H
#include <CGAL/Iso_rectangle_2.h>
#endif // CGAL_ISO_RECTANGLE_2_H
#ifndef CGAL_ISO_SQUARE_STATIC_2_H
#include <CGAL/Iso_square_static_2.h>
#endif // CGAL_ISO_SQUARE_STATIC_2_H

#line 46 "code_formatting.awi"
CGAL_BEGIN_NAMESPACE
#line 36 "pierce_traits.awi"

template < class _R >
struct Piercing_traits_cartesian {
  typedef typename _R::FT              FT;
  typedef CGAL::Point_2< _R >          Point_2;
  typedef CGAL::Iso_rectangle_2< _R >  Iso_rectangle_2;

  #line 65 "pierce_traits.awi"
  struct X : public CGAL_STD::unary_function< Point_2, FT >
  {
    FT
    operator()( const Point_2& p) const
    { return p.x(); }
  };
  
  struct Y : public CGAL_STD::unary_function< Point_2, FT >
  {
    FT
    operator()( const Point_2& p) const
    { return p.y(); }
  };
#line 44 "pierce_traits.awi"
  #line 81 "pierce_traits.awi"
  struct Xmin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmin(); }
  };
  
  struct Xmax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmax(); }
  };
  
  struct Ymin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymin(); }
  };
  
  struct Ymax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymax(); }
  };
#line 45 "pierce_traits.awi"
  #line 111 "pierce_traits.awi"
  struct Build_point
  : public CGAL_STD::binary_function< FT, FT, Point_2 >
  {
    Point_2
    operator()( const FT& px, const FT& py) const
    { return Point_2( px, py); }
  };
#line 46 "pierce_traits.awi"
  #line 138 "pierce_traits.awi"
  struct Build_rectangle
  : public CGAL_STD::unary_function< Point_2, Iso_rectangle_2 >
  {
    Iso_rectangle_2
    operator()( const Point_2& p) const
    { return Iso_rectangle_2( p, p); }
  };
}; // Piercing_traits_cartesian
#line 51 "pierce_traits.awi"
template < class _R >
struct Piercing_traits_homogeneous {
  typedef typename _R::FT              FT;
  typedef CGAL::Point_2< _R >          Point_2;
  typedef CGAL::Iso_rectangle_2< _R >  Iso_rectangle_2;

  #line 65 "pierce_traits.awi"
  struct X : public CGAL_STD::unary_function< Point_2, FT >
  {
    FT
    operator()( const Point_2& p) const
    { return p.x(); }
  };
  
  struct Y : public CGAL_STD::unary_function< Point_2, FT >
  {
    FT
    operator()( const Point_2& p) const
    { return p.y(); }
  };
#line 58 "pierce_traits.awi"
  #line 81 "pierce_traits.awi"
  struct Xmin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmin(); }
  };
  
  struct Xmax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmax(); }
  };
  
  struct Ymin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymin(); }
  };
  
  struct Ymax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymax(); }
  };
#line 59 "pierce_traits.awi"
  #line 121 "pierce_traits.awi"
  struct Build_point
  : public CGAL_STD::binary_function< FT, FT, Point_2 >
  {
    Point_2
    operator()( const FT& px, const FT& py) const
    {
      if ( px.denominator() == py.denominator())
        return Point_2( px.numerator(), py.numerator(), py.denominator());
      else
        return Point_2( px.numerator() * py.denominator(),
                        py.numerator() * px.denominator(),
                        px.denominator() * py.denominator());
    }
  };
#line 60 "pierce_traits.awi"
  #line 138 "pierce_traits.awi"
  struct Build_rectangle
  : public CGAL_STD::unary_function< Point_2, Iso_rectangle_2 >
  {
    Iso_rectangle_2
    operator()( const Point_2& p) const
    { return Iso_rectangle_2( p, p); }
  };
}; // Piercing_traits_homogeneous
#line 152 "pierce_traits.awi"
template < class _R >
struct Piercing_squares_traits_cartesian
: public Piercing_traits_cartesian< _R >
{
  typedef Piercing_traits_cartesian< _R >  Base;
  typedef Base::Point_2                    Point_2;
  typedef CGAL::Iso_square_static_2< _R >  Iso_rectangle_2;

  #line 81 "pierce_traits.awi"
  struct Xmin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmin(); }
  };
  
  struct Xmax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmax(); }
  };
  
  struct Ymin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymin(); }
  };
  
  struct Ymax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymax(); }
  };
#line 161 "pierce_traits.awi"
  #line 138 "pierce_traits.awi"
  struct Build_rectangle
  : public CGAL_STD::unary_function< Point_2, Iso_rectangle_2 >
  {
    Iso_rectangle_2
    operator()( const Point_2& p) const
    { return Iso_rectangle_2( p); }
  };
};
#line 166 "pierce_traits.awi"
template < class _R >
struct Piercing_squares_traits_homogeneous
: public Piercing_traits_homogeneous< _R >
{
  typedef Piercing_traits_homogeneous< _R >  Base;
  typedef Base::Point_2                      Point_2;
  typedef CGAL::Iso_square_static_2< _R >    Iso_rectangle_2;

  #line 81 "pierce_traits.awi"
  struct Xmin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmin(); }
  };
  
  struct Xmax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.xmax(); }
  };
  
  struct Ymin : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymin(); }
  };
  
  struct Ymax : public CGAL_STD::unary_function< Iso_rectangle_2, FT >
  {
    FT
    operator()( const Iso_rectangle_2& r) const
    { return r.ymax(); }
  };
#line 175 "pierce_traits.awi"
  #line 138 "pierce_traits.awi"
  struct Build_rectangle
  : public CGAL_STD::unary_function< Point_2, Iso_rectangle_2 >
  {
    Iso_rectangle_2
    operator()( const Point_2& p) const
    { return Iso_rectangle_2( p); }
  };
};
#line 256 "pierce_traits.awi"
#line 50 "code_formatting.awi"
CGAL_END_NAMESPACE
#line 257 "pierce_traits.awi"


#endif // ! (PIERCE_RECTANGLES_2_TRAITS_H)

#line 12 "code_formatting.awi"
// ----------------------------------------------------------------------------
// ** EOF
// ----------------------------------------------------------------------------

