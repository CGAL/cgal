// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// release       : 
// release_date  : 2000, August 16
//
// source        : webS2/S2.lw
// file          : include/CGAL/SimpleCartesian/Iso_rectangleS2.h
// package       : S2 (1.7)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 1.6
// revision_date : 27 Jun 2000
// author(s)     : Stefan Schirra
//                 based on code by
//                 Andreas Fabri and
//                 Herve Brönnimann
//
// coordinator   : MPI, Saarbrücken
// ======================================================================


#ifndef CGAL_ISO_RECTANGLES2_H
#define CGAL_ISO_RECTANGLES2_H

#include <CGAL/SimpleCartesian/PointS2.h>

CGAL_BEGIN_NAMESPACE

template <class FT>
class Iso_rectangleS2
{
public:
                  Iso_rectangleS2() {}
                  Iso_rectangleS2(const PointS2<FT>& p, const PointS2<FT>& q);

  bool            operator==(const Iso_rectangleS2<FT>& s) const;
  bool            operator!=(const Iso_rectangleS2<FT>& s) const;
  int             id() const;

  const PointS2<FT>&     min() const;
  const PointS2<FT>&     max() const;
  PointS2<FT>    vertex(int i) const;
  PointS2<FT>    operator[](int i) const;

  Iso_rectangleS2<FT> transform(const Aff_transformationS2<FT>& t) const;

  Bounded_side    bounded_side(const PointS2<FT>& p) const;

  bool            has_on_boundary(const PointS2<FT>& p) const;

  bool            has_on_bounded_side(const PointS2<FT>& p) const;
  bool            has_on_unbounded_side(const PointS2<FT>& p) const;

  bool            is_degenerate() const;

  Bbox_2          bbox() const;

  FT              xmin() const;
  FT              ymin() const;
  FT              xmax() const;
  FT              ymax() const;

// private:
  PointS2<FT>  e0;
  PointS2<FT>  e1;
};


template < class FT >
inline
Iso_rectangleS2<FT>::Iso_rectangleS2(const PointS2<FT>& p,
                                       const PointS2<FT>& q)
{
  FT vx0 = p.x();
  FT vy0 = p.y();
  FT vx1 = q.x();
  FT vy1 = q.y();

  bool b1 = false,
       b2 = false;
  if ( (b1 = (vx0 > vx1)) || (b2 = (vy0 > vy1)) ) 
  {
    if (b1 && b2) 
    { e0 = q; e1 = p; }
    else 
    {
      if (vx0 > vx1)
      {
        FT z = vx1;
        vx1 = vx0;
        vx0 = z;
      }
      if (vy0 > vy1)
      {
        FT z = vy1;
        vy1 = vy0;
        vy0 = z;
      }
      e0 = PointS2<FT>(vx0,vy0); 
      e1 = PointS2<FT>(vx1,vy1);
    }
  }
  else 
  { e0 = p; e1 = q; }
}


template < class FT >
inline
bool 
Iso_rectangleS2<FT>::operator==(const Iso_rectangleS2<FT>& r) const
{ return  vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2); }

template < class FT >
inline
bool 
Iso_rectangleS2<FT>::operator!=(const Iso_rectangleS2<FT>& r) const
{ return !(*this == r); }

template < class FT >
inline
const PointS2<FT>&
Iso_rectangleS2<FT>::min() const
{ return  e0; }

template < class FT >
inline
const PointS2<FT>& 
Iso_rectangleS2<FT>::max() const
{ return  e1; }

template < class FT >
inline 
FT 
Iso_rectangleS2<FT>::xmin() const
{ return  min().x(); }

template < class FT >
inline 
FT 
Iso_rectangleS2<FT>::ymin() const
{ return  min().y(); }

template < class FT >
inline 
FT 
Iso_rectangleS2<FT>::xmax() const
{ return  max().x(); }

template < class FT >
inline 
FT 
Iso_rectangleS2<FT>::ymax() const
{ return  max().y(); }

template < class FT >
PointS2<FT>
Iso_rectangleS2<FT>::vertex(int i) const
{
  switch (i%4) {
  case 0: return min();
  case 1: return PointS2<FT>(xmax(), ymin());
  case 2: return max();
  default: return PointS2<FT>(xmin(), ymax());
  }
}

template < class FT >
inline 
PointS2<FT>
Iso_rectangleS2<FT>::operator[](int i) const
{ return vertex(i); }

template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side Iso_rectangleS2<FT>::bounded_side(const PointS2<FT>& p) const
{
  bool x_incr = (xmin() < p.x()) &&  (p.x() < xmax()),
       y_incr = (ymin() < p.y()) &&  (p.y() < ymax());
  if( x_incr )
    {
      if( y_incr )
        {
          return ON_BOUNDED_SIDE;
        }
      if( (p.y() == ymin()) || (ymax() == p.y()) )
        {
          return ON_BOUNDARY;
        }
    }
  if( (p.x() == xmin()) || (xmax() == p.x()) )
    {
      if( y_incr || (p.y() == ymin()) || (ymax() == p.y()) )
        {
          return ON_BOUNDARY;
        }
    }

  return ON_UNBOUNDED_SIDE;
}

template < class FT >
inline
bool 
Iso_rectangleS2<FT>::has_on_boundary(const PointS2<FT>& p) const
{ return bounded_side(p) == ON_BOUNDARY; }

template < class FT >
inline bool Iso_rectangleS2<FT>::has_on_bounded_side(
                                   const PointS2<FT>& p) const
{ return bounded_side(p) == ON_BOUNDED_SIDE; }

template < class FT >
inline bool Iso_rectangleS2<FT>::has_on_unbounded_side(
                                   const PointS2<FT>& p) const
{ return bounded_side(p) == ON_UNBOUNDED_SIDE; }

template < class FT >
bool Iso_rectangleS2<FT>::is_degenerate() const
{ return (xmin() == xmax()) || (ymin() ==ymax()); }

template < class FT >
inline
Bbox_2 Iso_rectangleS2<FT>::bbox() const
{
  return Bbox_2(CGAL::to_double(xmin()), CGAL::to_double(ymin()),
                CGAL::to_double(xmax()), CGAL::to_double(ymax()));
}

template < class FT >
inline
Iso_rectangleS2<FT> 
Iso_rectangleS2<FT>::transform(const Aff_transformationS2<FT>& t) const
{
  return Iso_rectangleS2<FT>(t.transform(vertex(0)),
                             t.transform(vertex(2)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLES2
template < class FT >
std::ostream& 
operator<<(std::ostream& os, const Iso_rectangleS2<FT> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r[0] << ' ' << r[2];
    case IO::BINARY :
        return os << r[0] << r[2];
    default:
        return os << "Iso_rectangleS2(" << r[0] << ", " << r[2] << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLES2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLES2
template < class FT >
CGAL_KERNEL_MEDIUM_INLINE
std::istream& 
operator>>(std::istream& is, Iso_rectangleS2<FT> &r)
{
    PointS2<FT> p, q;

    is >> p >> q;

    r = Iso_rectangleS2<FT>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLES2



CGAL_END_NAMESPACE

#endif
