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
// file          : Iso_cuboidH3.h
// package       : H3
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Stefan Schirra
//
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_ISO_CUBOIDH3_H
#define CGAL_ISO_CUBOIDH3_H

#include <CGAL/Twotuple.h>
#include <CGAL/PointH3.h>
#include <CGAL/predicates_on_pointsH3.h>

CGAL_BEGIN_NAMESPACE

template <class R_>
class Iso_cuboidH3
  : public R_::Iso_cuboid_handle_3
{
public:
  typedef R_                R;
  typedef typename R::RT    RT;
  typedef typename R::FT    FT;

  typedef typename R::Iso_cuboid_handle_3  Iso_cuboid_handle_3_;
  typedef typename Iso_cuboid_handle_3_::element_type Iso_cuboid_ref_3;

  Iso_cuboidH3()
    : Iso_cuboid_handle_3_(Iso_cuboid_ref_3()) {}

  Iso_cuboidH3(const PointH3<R>& p, const PointH3<R>& q);

  Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz, 
               const RT& hw);

  Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
               const RT& max_hx, const RT& max_hy, const RT& max_hz);

  bool      operator==(const Iso_cuboidH3<R>& s) const;
  bool      operator!=(const Iso_cuboidH3<R>& s) const;

  PointH3<R>  min() const;
  PointH3<R>  max() const;
  PointH3<R>  vertex(int i) const;
  PointH3<R>  operator[](int i) const;

  Iso_cuboidH3<R>
            transform(const Aff_transformationH3<R>& t) const;
  Bounded_side
            bounded_side(const PointH3<R>& p) const;
  bool      has_on(const PointH3<R>& p) const;
  bool      has_on_boundary(const PointH3<R>& p) const;
  bool      has_on_bounded_side(const PointH3<R>& p) const;
  bool      has_on_unbounded_side(const PointH3<R>& p) const;
  bool      is_degenerate() const;
  Bbox_2    bbox() const;
  FT        xmin() const;
  FT        ymin() const;
  FT        zmin() const;
  FT        xmax() const;
  FT        ymax() const;
  FT        zmax() const;
  FT        min_coord(int i) const;
  FT        max_coord(int i) const;

  FT        volume() const;

};

template < class R >
CGAL_KERNEL_CTOR_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const PointH3<R>& p, const PointH3<R>& q)
{
  bool px_k_qx = ( p.hx()*q.hw() < q.hx()*p.hw() );
  bool py_k_qy = ( p.hy()*q.hw() < q.hy()*p.hw() );
  bool pz_k_qz = ( p.hz()*q.hw() < q.hz()*p.hw() );

  RT minx;
  RT miny;
  RT minz;
  RT maxx;
  RT maxy;
  RT maxz;
  RT minw = p.hw()*q.hw();
  RT maxw = p.hw()*q.hw();
  if ( px_k_qx )
  {
      minx = p.hx()*q.hw();
      maxx = q.hx()*p.hw();
  }
  else
  {
      minx = q.hx()*p.hw();
      maxx = p.hx()*q.hw();
  }
  if ( py_k_qy )
  {
      miny = p.hy()*q.hw();
      maxy = q.hy()*p.hw();
  }
  else
  {
      miny = q.hy()*p.hw();
      maxy = p.hy()*q.hw();
  }
  if ( pz_k_qz )
  {
      minz = p.hz()*q.hw();
      maxz = q.hz()*p.hw();
  }
  else
  {
      minz = q.hz()*p.hw();
      maxz = p.hz()*q.hw();
  }
  initialize_with( Iso_cuboid_ref_3 ( PointH3<R>(minx, miny, minz, minw),
                                      PointH3<R>(maxx, maxy, maxz, maxw) ));
}

template < class R >
CGAL_KERNEL_CTOR_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
             const RT& max_hx, const RT& max_hy, const RT& max_hz)
{
  initialize_with( 
     Iso_cuboid_ref_3 ( PointH3<R>(min_hx, min_hy, min_hz, RT(1)),
                        PointH3<R>(max_hx, max_hy, max_hz, RT(1)) ));
}

template < class R >
CGAL_KERNEL_CTOR_LARGE_INLINE
Iso_cuboidH3<R>::
Iso_cuboidH3(const RT& min_hx, const RT& min_hy, const RT& min_hz,
             const RT& max_hx, const RT& max_hy, const RT& max_hz, 
             const RT& hw)
{
  initialize_with( Iso_cuboid_ref_3( PointH3<R>(min_hx, min_hy, min_hz, hw),
                                     PointH3<R>(max_hx, max_hy, max_hz, hw) ));
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::
operator==(const Iso_cuboidH3<R>& r) const
{ return  (min() == r.min()) && (max() == r.max()); }

template < class R >
inline
bool
Iso_cuboidH3<R>::
operator!=(const Iso_cuboidH3<R>& r) const
{ return !(*this == r); }

template < class R >
inline
PointH3<R>
Iso_cuboidH3<R>::min() const
{ return  Ptr()->e0; }

template < class R >
inline
PointH3<R>
Iso_cuboidH3<R>::max() const
{ return  Ptr()->e1; }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::xmin() const
{ return  FT( min().hx() ) / FT( min().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::ymin() const
{ return  FT( min().hy() ) / FT( min().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::zmin() const
{ return  FT( min().hz() ) / FT( min().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::xmax() const
{ return  FT( max().hx() ) / FT( max().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::ymax() const
{ return  FT( max().hy() ) / FT( max().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::zmax() const
{ return  FT( max().hz() ) / FT( max().hw() ); }

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::min_coord(int i) const
{ 
   CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
   if ( i == 0 )
       return xmin();
   else if (i == 1)
       return ymin();
   return zmin();
}

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::max_coord(int i) const
{ 
   CGAL_kernel_precondition(i == 0 || i == 1 || i == 2);
   if ( i == 0 )
      return xmax();
   else if ( i == 1 )
      return ymax();
   return zmax();
}

template < class R >
inline
typename Iso_cuboidH3<R>::FT
Iso_cuboidH3<R>::volume() const
{ return  (xmax() - xmin()) * (ymax() - ymin()) * (zmax() - zmin()); }

template < class R >
CGAL_KERNEL_LARGE_INLINE
PointH3<R>
Iso_cuboidH3<R>::vertex(int i) const
{
  switch (i%8)
  {
    case 0: return min();
    case 1: return PointH3<R>( max().hx(), min().hy(),
		                   min().hz(), min().hw() );
    case 2: return PointH3<R>( max().hx(), max().hy(),
		                   min().hz(), min().hw() );
    case 3: return PointH3<R>( min().hx(), max().hy(),
		                   min().hz(), min().hw() );
    case 4: return PointH3<R>( min().hx(), max().hy(),
		                   max().hz(), min().hw() );
    case 5: return PointH3<R>( min().hx(), min().hy(),
		                   max().hz(), min().hw() );
    case 6: return PointH3<R>( max().hx(), min().hy(),
		                   max().hz(), min().hw() );
    case 7: return max();
  }
  return PointH3<R>();
}

template < class R >
inline
PointH3<R>
Iso_cuboidH3<R>::operator[](int i) const
{ return vertex(i); }

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidH3<R>::
bounded_side(const PointH3<R>& p) const
{
  if (    (p.hx()*min().hw() < min().hx()*p.hw() )
        ||(p.hy()*min().hw() < min().hy()*p.hw() )
        ||(p.hz()*min().hw() < min().hz()*p.hw() )
        ||(p.hx()*max().hw() > max().hx()*p.hw() )
        ||(p.hy()*max().hw() > max().hy()*p.hw() )
        ||(p.hz()*max().hw() > max().hz()*p.hw() )  )
  { return ON_UNBOUNDED_SIDE; }
  if (    (p.hx()*min().hw() == min().hx()*p.hw() )
        ||(p.hy()*min().hw() == min().hy()*p.hw() )
        ||(p.hz()*min().hw() == min().hz()*p.hw() )
        ||(p.hx()*max().hw() == max().hx()*p.hw() )
        ||(p.hy()*max().hw() == max().hy()*p.hw() )
        ||(p.hz()*max().hw() == max().hz()*p.hw() )  )
  { return ON_BOUNDARY; }
  else
  { return ON_BOUNDED_SIDE; }
}

template < class R >
inline
bool
Iso_cuboidH3<R>::has_on_boundary(const PointH3<R>& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_cuboidH3<R>::has_on(const PointH3<R>& p) const
{ return ( bounded_side(p) == ON_BOUNDARY ); }

template < class R >
inline
bool
Iso_cuboidH3<R>::
has_on_bounded_side(const PointH3<R>& p) const
{ return ( bounded_side(p) == ON_BOUNDED_SIDE ); }

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::
has_on_unbounded_side(const PointH3<R>& p) const
{
  return (   ( lexicographically_xyz_smaller(p,min() ))
           ||( lexicographically_xyz_smaller(max(),p ))  );
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidH3<R>::is_degenerate() const
{
  return (  ( min().hx() == max().hx() )
         || ( min().hy() == max().hy() )
         || ( min().hz() == max().hz() ) );
}

template < class R >
inline
Bbox_2
Iso_cuboidH3<R>::bbox() const
{ return  min().bbox() + max().bbox(); }

template < class R >
CGAL_KERNEL_INLINE
Iso_cuboidH3<R>
Iso_cuboidH3<R>::
transform(const Aff_transformationH3<R>&t) const
{
  return Iso_cuboidH3<R>(t.transform(min() ),
                             t.transform(max() ) );
}

#ifndef CGAL_NO_OSTREAM_INSERT_ISO_CUBOIDH3
template < class R >
std::ostream& operator<<(std::ostream& os, const Iso_cuboidH3<R>& r)
{
  switch(os.iword(IO::mode))
  {
    case IO::ASCII :
        return os << r.min() << ' ' << r.max();
    case IO::BINARY :
        return os << r.min() << r.max();
    default:
        return os << "Iso_cuboidH3(" << r.min() << ", " << r.max() << ")";
  }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_CUBOIDH3

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOIDH3
template < class R >
std::istream& operator>>(std::istream& is, Iso_cuboidH3<R>& r)
{
  PointH3<R> p, q;
  is >> p >> q;
  r = Iso_cuboidH3<R>(p, q);
  return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_CUBOIDH3

CGAL_END_NAMESPACE

#endif // CGAL_ISO_CUBOIDH3_H
