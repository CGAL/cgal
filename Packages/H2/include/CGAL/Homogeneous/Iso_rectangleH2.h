// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 
// 
// file          : include/CGAL/Homogeneous/Iso_rectangleH2.h
// package       : H2
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 
#ifndef CGAL_ISO_RECTANGLEH2_H
#define CGAL_ISO_RECTANGLEH2_H

#include <CGAL/Twotuple.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_rectangleH2
  : public R_::template Handle<Twotuple<typename R_::Point_2> >::type
{
CGAL_VC7_BUG_PROTECTED
  typedef typename R_::FT                   FT;
  typedef typename R_::RT                   RT;
  typedef typename R_::Point_2              Point_2;
  typedef typename R_::Aff_transformation_2 Aff_transformation_2;

  typedef Twotuple<Point_2>                        rep;
  typedef typename R_::template Handle<rep>::type  base;

public:
  typedef R_                                    R;

  Iso_rectangleH2()
    : base(rep()) {}

  Iso_rectangleH2(const Point_2& p, const Point_2& q);

  Iso_rectangleH2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy);

  Iso_rectangleH2(const RT& min_hx, const RT& min_hy, 
                  const RT& max_hx, const RT& max_hy, const RT& hw);

  bool      operator==(const Iso_rectangleH2<R>& s) const;
  bool      operator!=(const Iso_rectangleH2<R>& s) const;

  const Point_2 & min() const;
  const Point_2 & max() const;
  Point_2  vertex(int i) const;
  Point_2  operator[](int i) const;

  Iso_rectangleH2<R>
            transform(const Aff_transformation_2& t) const;

  Bounded_side
            bounded_side(const Point_2& p) const;
  bool      has_on(const Point_2& p) const;
  bool      has_on_boundary(const Point_2& p) const;
  bool      has_on_bounded_side(const Point_2& p) const;
  bool      has_on_unbounded_side(const Point_2& p) const;
  bool      is_degenerate() const;

  Bbox_2    bbox() const;

  FT        xmin() const;
  FT        ymin() const;
  FT        xmax() const;
  FT        ymax() const;
  FT        min_coord(int i) const;
  FT        max_coord(int i) const;

  FT        area() const;
};

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Iso_rectangleH2<R>::
Iso_rectangleH2(const typename Iso_rectangleH2<R>::Point_2& p,
	        const typename Iso_rectangleH2<R>::Point_2& q)
{
  bool px_g_qx = ( p.hx()*q.hw() > q.hx()*p.hw() );
  bool py_g_qy = ( p.hy()*q.hw() > q.hy()*p.hw() );
  if ( px_g_qx || py_g_qy)
  {
      if ( px_g_qx && py_g_qy )
      {
          initialize_with( rep(q,p) );
      }
      else
      {
         if ( px_g_qx )
         {
             initialize_with( rep(
             Point_2(q.hx()*p.hw(), p.hy()*q.hw(), q.hw()*p.hw() ),
             Point_2(p.hx()*q.hw(), q.hy()*p.hw(), q.hw()*p.hw() )) );
         }
         if ( py_g_qy )
         {
             initialize_with( rep(
             Point_2(p.hx()*q.hw(), q.hy()*p.hw(), q.hw()*p.hw() ),
             Point_2(q.hx()*p.hw(), p.hy()*q.hw(), q.hw()*p.hw() )) );
         }
      }
  }
  else
  {
      initialize_with( rep(p,q) );
  }
}

template < class R >
inline
Iso_rectangleH2<R>::Iso_rectangleH2(const RT& min_hx, const RT& min_hy, 
                                    const RT& max_hx, const RT& max_hy)
{
  initialize_with( rep( Point_2(min_hx, min_hy), 
                        Point_2(max_hx, max_hy)) );
}

template < class R >
inline
Iso_rectangleH2<R>::Iso_rectangleH2(const RT& min_hx, const RT& min_hy, 
                                    const RT& max_hx, const RT& max_hy,
                                    const RT& hw)
{
  initialize_with( rep( Point_2(min_hx, min_hy, hw), 
                        Point_2(max_hx, max_hy, hw)) );
}

template < class R >
inline
bool
Iso_rectangleH2<R>::operator==(const Iso_rectangleH2<R>& r) const
{ return  vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2); }

template < class R >
inline
bool
Iso_rectangleH2<R>::operator!=(const Iso_rectangleH2<R>& r) const
{ return !(*this == r); }

template < class R >
inline
const typename Iso_rectangleH2<R>::Point_2 &
Iso_rectangleH2<R>::min() const
{ return  Ptr()->e0; }

template < class R >
inline
const typename Iso_rectangleH2<R>::Point_2 &
Iso_rectangleH2<R>::max() const
{ return  Ptr()->e1; }

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::xmin() const
{ return  FT( min().hx() ) / FT( min().hw() ); }

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::ymin() const
{ return  FT( min().hy() ) / FT( min().hw() ); }

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::xmax() const
{ return  FT( max().hx() ) / FT( max().hw() ); }

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::ymax() const
{ return  FT( max().hy() ) / FT( max().hw() ); }

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::min_coord(int i) const
{ 
   CGAL_kernel_precondition ( i == 0 || i == 1 );
   if (i == 0)
      return xmin();
   else
      return ymin();
}

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::max_coord(int i) const
{ 
   CGAL_kernel_precondition ( i == 0 || i == 1 );
   if (i == 0)
      return xmax();
   else
      return ymax();
}

template < class R >
inline
typename Iso_rectangleH2<R>::FT
Iso_rectangleH2<R>::area() const
{ return  (xmax() - xmin()) * (ymax() - ymin()); }

template < class R >
CGAL_KERNEL_INLINE
typename Iso_rectangleH2<R>::Point_2
Iso_rectangleH2<R>::vertex(int i) const
{
  switch (i%4)
  {
    case 0:
        return min();
    case 1:
        return Point_2( max().hx()*min().hw(),
                           min().hy()*max().hw(),
                           min().hw()*max().hw() );
    case 2:
        return max();
    default: // case 3:
        return Point_2( min().hx()*max().hw(),
                           max().hy()*min().hw(),
                           min().hw()*max().hw() );
  }
}

template < class R >
inline
typename Iso_rectangleH2<R>::Point_2
Iso_rectangleH2<R>::operator[](int i) const
{ return vertex(i); }

template < class R >
CGAL_KERNEL_INLINE
Bounded_side
Iso_rectangleH2<R>::
bounded_side(const typename Iso_rectangleH2<R>::Point_2& p) const
{
  Oriented_side wrt_min = _where_wrt_L_wedge(min(),p);
  Oriented_side wrt_max = _where_wrt_L_wedge(p,max());
  if (( wrt_min == ON_NEGATIVE_SIDE )||( wrt_max == ON_NEGATIVE_SIDE))
  {
      return ON_UNBOUNDED_SIDE;
  }
  if (  ( wrt_min == ON_ORIENTED_BOUNDARY )
      ||( wrt_max == ON_ORIENTED_BOUNDARY ) )
  {
      return ON_BOUNDARY;
  }
  return ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleH2<R>::
has_on_boundary(const typename Iso_rectangleH2<R>::Point_2& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_rectangleH2<R>::has_on(const typename Iso_rectangleH2<R>::Point_2& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_rectangleH2<R>::
has_on_bounded_side(const typename Iso_rectangleH2<R>::Point_2& p) const
{ return ( bounded_side(p) == ON_BOUNDED_SIDE ); }

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_rectangleH2<R>::
has_on_unbounded_side(const typename Iso_rectangleH2<R>::Point_2& p) const
{
  return (  (_where_wrt_L_wedge(min(),p) == ON_NEGATIVE_SIDE)
          ||(_where_wrt_L_wedge(p,max()) == ON_NEGATIVE_SIDE) );
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_rectangleH2<R>::is_degenerate() const
{
 return (   ( min().hx()*max().hw() == max().hx()*min().hw() )
         || ( min().hy()*max().hw() == max().hy()*min().hw() ) );
}
template < class R >
inline
Bbox_2
Iso_rectangleH2<R>::bbox() const
{ return  min().bbox() + max().bbox(); }

template < class R >
CGAL_KERNEL_INLINE
Iso_rectangleH2<R>
Iso_rectangleH2<R>::
transform(const typename Iso_rectangleH2<R>::Aff_transformation_2&t) const
{
  return Iso_rectangleH2<R>(t.transform(min() ), t.transform(max() ) );
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEH2
template < class R >
std::ostream& operator<<(std::ostream& os, const Iso_rectangleH2<R>& r)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << r[0] << ' ' << r[2];
    case IO::BINARY :
        return os << r[0] << r[2];
    default:
        return os << "Iso_rectangleH2(" << r[0] << ", " << r[2] << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEH2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEH2
template < class R >
std::istream& operator>>(std::istream& is, Iso_rectangleH2<R>& r)
{
  typename R::Point_2 p, q;
  is >> p >> q;
  r = Iso_rectangleH2<R>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEH2

CGAL_END_NAMESPACE

#endif // CGAL_ISO_RECTANGLEH2_H
