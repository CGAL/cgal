// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann

#ifndef CGAL_CARTESIAN_ISO_CUBOID_3_C
#define CGAL_CARTESIAN_ISO_CUBOID_3_C

#include <CGAL/Bbox_3.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>

#ifndef CGAL_CARTESIANR_EDEFINE_NAMES_3_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
inline
_Twotuple< typename Iso_cuboidC3<R CGAL_CTAG>::Point_3 > *
Iso_cuboidC3<R CGAL_CTAG>::ptr() const
{
  return (_Twotuple<Point_3>*)PTR;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
Iso_cuboidC3<R CGAL_CTAG>::Iso_cuboidC3()
{
  PTR = new _Twotuple< Point_3 >;
}

template < class R >
CGAL_KERNEL_CTOR_INLINE
Iso_cuboidC3<R CGAL_CTAG>::Iso_cuboidC3(const Iso_cuboidC3<R CGAL_CTAG>& r)
  : Handle((Handle&)r)
{}

template < class R >
CGAL_KERNEL_CTOR_LARGE_INLINE
Iso_cuboidC3<R CGAL_CTAG>::
Iso_cuboidC3(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p,
             const Iso_cuboidC3<R CGAL_CTAG>::Point_3& q)
{
  FT minx, maxx, miny, maxy, minz, maxz;
  if (p.x() < q.x()) { minx = p.x(); maxx = q.x(); }
  else               { minx = q.x(); maxx = p.x(); }
  if (p.y() < q.y()) { miny = p.y(); maxy = q.y(); }
  else               { miny = q.y(); maxy = p.y(); }
  if (p.z() < q.z()) { minz = p.z(); maxz = q.z(); }
  else               { minz = q.z(); maxz = p.z(); }
  PTR = new _Twotuple<Point_3 >(Point_3(minx, miny, minz),
                                Point_3(maxx, maxy, maxz) );
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::~Iso_cuboidC3()
{}

template < class R >
CGAL_KERNEL_INLINE
Iso_cuboidC3<R CGAL_CTAG>&
Iso_cuboidC3<R CGAL_CTAG>::operator=(const Iso_cuboidC3<R CGAL_CTAG>& r)
{
 Handle::operator=(r);
 return *this;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::operator==(const Iso_cuboidC3<R CGAL_CTAG>& r) const
{
  return min() == r.min() && max() == r.max();
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::operator!=(const Iso_cuboidC3<R CGAL_CTAG>& r) const
{
  return !(*this == r);
}

template < class R >
inline
int
Iso_cuboidC3<R CGAL_CTAG>::id() const
{
  return (int)PTR;
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::min() const
{
  return  ptr()->e0;
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::max() const
{
  return  ptr()->e1;
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::xmin() const
{
  return min().x();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::ymin() const
{
  return min().y();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::zmin() const
{
  return min().z();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::xmax() const
{
  return max().x();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::ymax() const
{
  return max().y();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::FT
Iso_cuboidC3<R CGAL_CTAG>::zmax() const
{
  return max().z();
}

template < class R >
CGAL_KERNEL_LARGE_INLINE
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::vertex(int i) const
{
  switch (i%8)
  {
    case 0: return min();
    case 1: return Point_3(max().hx(), min().hy(), min().hz());
    case 2: return Point_3(max().hx(), max().hy(), min().hz());
    case 3: return Point_3(min().hx(), max().hy(), min().hz());
    case 4: return Point_3(min().hx(), max().hy(), max().hz());
    case 5: return Point_3(min().hx(), min().hy(), max().hz());
    case 6: return Point_3(max().hx(), min().hy(), max().hz());
    case 7: return max();
  }
  return Point_3();
}

template < class R >
inline
Iso_cuboidC3<R CGAL_CTAG>::Point_3
Iso_cuboidC3<R CGAL_CTAG>::operator[](int i) const
{
  return vertex(i);
}

template < class R >
CGAL_KERNEL_MEDIUM_INLINE
Bounded_side
Iso_cuboidC3<R CGAL_CTAG>::
bounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  if (strict_dominance(p,min()) && strict_dominance(max(),p) )
    return ON_BOUNDED_SIDE;
  if (dominance(p,min()) && dominance(max(),p))
    return ON_BOUNDARY;
  return ON_UNBOUNDED_SIDE;
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on_boundary(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDARY;
}

template < class R >
inline
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on_bounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return bounded_side(p) == ON_BOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::
has_on_unbounded_side(const Iso_cuboidC3<R CGAL_CTAG>::Point_3& p) const
{
  return  bounded_side(p) == ON_UNBOUNDED_SIDE;
}

template < class R >
CGAL_KERNEL_INLINE
bool
Iso_cuboidC3<R CGAL_CTAG>::is_degenerate() const
{
  return min().hx() == max().hx()
      || min().hy() == max().hy()
      || min().hz() == max().hz();
}

template < class R >
inline
Bbox_3
Iso_cuboidC3<R CGAL_CTAG>::bbox() const
{
  return min().bbox() + max().bbox();
}

template < class R >
CGAL_KERNEL_INLINE
Iso_cuboidC3<R CGAL_CTAG>
Iso_cuboidC3<R CGAL_CTAG>::
transform(const Iso_cuboidC3<R CGAL_CTAG>::Aff_transformation_3&t) const
{
  return Self(t.transform(min()), t.transform(max()) );
}

#ifndef NO_OSTREAM_INSERT_ISO_CUBOIDC3
template < class R >
std::ostream &
operator<<(std::ostream& os, const Iso_cuboidC3<R CGAL_CTAG>& r)
{
  switch(os.iword(IO::mode)) {
  case IO::ASCII :
    return os << min() << ' ' << max();
  case IO::BINARY :
    return os << min() << max();
  default:
    return os << "Iso_cuboidC3(" << min() << ", " << max() << ")";
  }
}
#endif // NO_OSTREAM_INSERT_ISO_CUBOIDC3

#ifndef NO_ISTREAM_EXTRACT_ISO_CUBOIDC3
template < class R >
std::istream &
operator>>(std::istream& is, Iso_cuboidC3<R CGAL_CTAG>& r)
{
  Iso_cuboidC3<R CGAL_CTAG>::Point_3 p, q;
  is >> p >> q;
  r = Iso_cuboidC3<R CGAL_CTAG>(p, q);
  return is;
}
#endif // NO_ISTREAM_EXTRACT_ISO_CUBOIDC3

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_ISO_CUBOID_3_C
