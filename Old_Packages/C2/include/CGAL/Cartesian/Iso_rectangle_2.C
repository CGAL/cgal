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
// release       :
// release_date  :
//
// file          : include/CGAL/Cartesian/Iso_rectangle_2.h
// source        : include/CGAL/Cartesian/Cartesian_2/Iso_rectangle_2.h
// revision      : include/CGAL/Cartesian/Cartesian_2/Iso_rectangle_2.h
// revision_date : $Date$
// author(s)     : Andreas.Fabri@sophia.inria.fr
//                 Herve Bronnimann@sophia.inria.fr
//
// coordinator   : INRIA Sophia-Antipolis (Herve.Bronnimann@sophia.inria.fr)
//
// ============================================================================

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef _MSC_VER
#define typename
#endif

#ifndef CGAL_CARTESIAN_ISO_RECTANGLE_2_C
#define CGAL_CARTESIAN_ISO_RECTANGLE_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
inline _Twotuple< typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 > *
Iso_rectangleC2<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<Point_2>*)PTR;
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::Iso_rectangleC2()
{
  PTR = new _Twotuple<Point_2>;
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::
Iso_rectangleC2(const Iso_rectangleC2<R CGAL_CTAG> &r)
  : Handle((Handle&)r)
{
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::
Iso_rectangleC2(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p,
                const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &q)
{
  FT vx0 = p.x();
  FT vy0 = p.y();
  FT vx1 = q.x();
  FT vy1 = q.y();

  bool b1 = false,
       b2 = false;
  if( (b1 = (vx0 > vx1)) || (b2 = (vy0 > vy1)) ) {
    if (b1 && b2) {
      PTR = new _Twotuple<Point_2>(q,p);
    } else {
      if (vx0 > vx1){
        FT z = vx1;
        vx1 = vx0;
        vx0 = z;
      }
      if (vy0 > vy1){
        FT z = vy1;
        vy1 = vy0;
        vy0 = z;
      }

      PTR = new _Twotuple<Point_2>(Point_2(vx0,vy0), Point_2(vx1,vy1));
    }
  }else {
    PTR = new _Twotuple<Point_2>(p,q);
  }
}


template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>::~Iso_rectangleC2()
{}


template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG> &
Iso_rectangleC2<R CGAL_CTAG>::operator=(const Iso_rectangleC2<R CGAL_CTAG> &r)
{

  Handle::operator=(r);
  return *this;
}
template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
operator==(const Iso_rectangleC2<R CGAL_CTAG> &r) const
{
  return  vertex(0) == r.vertex(0) && vertex(2) == r.vertex(2);
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
operator!=(const Iso_rectangleC2<R CGAL_CTAG> &r) const
{
  return !(*this == r);
}

template < class R >
inline
int Iso_rectangleC2<R CGAL_CTAG>::id() const
{
  return (int)PTR;
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::min() const
{
  return  ptr()->e0;
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::max() const
{
  return  ptr()->e1;
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT Iso_rectangleC2<R CGAL_CTAG>::xmin() const
{
  return  min().x();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT Iso_rectangleC2<R CGAL_CTAG>::ymin() const
{
  return  min().y();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT Iso_rectangleC2<R CGAL_CTAG>::xmax() const
{
  return  max().x();
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::FT Iso_rectangleC2<R CGAL_CTAG>::ymax() const
{
  return  max().y();
}



template < class R >
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::vertex(int i) const
{
  switch (i%4) {
  case 0: return min();
  case 1: return Point_2(xmax(), ymin());
  case 2: return max();
  default: return Point_2(xmin(), ymax());
  }
}

template < class R >
inline
typename Iso_rectangleC2<R CGAL_CTAG>::Point_2
Iso_rectangleC2<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_rectangleC2<R CGAL_CTAG>::
bounded_side(const Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
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

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_boundary(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_bounded_side(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_rectangleC2<R CGAL_CTAG>::
has_on_unbounded_side(const typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 &p) const
{
  return bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
bool Iso_rectangleC2<R CGAL_CTAG>::is_degenerate() const
{
  return (xmin() == xmax()) || (ymin() ==ymax());
}
template < class R >
inline
Bbox_2 Iso_rectangleC2<R CGAL_CTAG>::bbox() const
{
  return Bbox_2(CGAL::to_double(xmin()), CGAL::to_double(ymin()),
                CGAL::to_double(xmax()), CGAL::to_double(ymax()));
}

template < class R >
inline
Iso_rectangleC2<R CGAL_CTAG>
Iso_rectangleC2<R CGAL_CTAG>::
transform(const typename Iso_rectangleC2<R CGAL_CTAG>::Aff_transformation_2 &t) const
{
  // We need a precondition like this!!!
  // CGAL_kernel_precondition(t.is_axis_preserving());
  return Iso_rectangleC2<R CGAL_CTAG>(t.transform(vertex(0)),
                             t.transform(vertex(2)));
}


#ifndef CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEC2
template < class R >
std::ostream &
operator<<(std::ostream &os, const Iso_rectangleC2<R CGAL_CTAG> &r)
{
    switch(os.iword(IO::mode)) {
    case IO::ASCII :
        return os << r[0] << ' ' << r[2];
    case IO::BINARY :
        return os << r[0] << r[2];
    default:
        return os << "Iso_rectangleC2(" << r[0] << ", " << r[2] << ")";
    }
}
#endif // CGAL_NO_OSTREAM_INSERT_ISO_RECTANGLEC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEC2
template < class R >
CGAL_KERNEL_MEDIUM_INLINE
std::istream &
operator>>(std::istream &is, Iso_rectangleC2<R CGAL_CTAG> &r)
{
    typename Iso_rectangleC2<R CGAL_CTAG>::Point_2 p, q;

    is >> p >> q;

    r = Iso_rectangleC2<R CGAL_CTAG>(p, q);
    return is;
}
#endif // CGAL_NO_ISTREAM_EXTRACT_ISO_RECTANGLEC2

CGAL_END_NAMESPACE

#ifdef _MSC_VER
#undef typename
#endif

#endif
