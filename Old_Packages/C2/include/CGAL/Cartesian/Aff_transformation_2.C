#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C

#ifndef CGAL_DETERMINANT_H
#include <CGAL/determinant.h>
#endif // CGAL_DETERMINANT_H

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_REP_2_C
#include <CGAL/Cartesian/Aff_transformation_rep_2.C>
#endif // CGAL_CARTESIAN_+AFF_TRANSFORMATION_REP_2_C

CGAL_BEGIN_NAMESPACE

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2()
{
  PTR = new Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1));
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Identity)
{
  PTR = new Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1));
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::Aff_transformationC2
  (const Aff_transformationC2<R CGAL_CTAG> &t)
  : Handle(t)
{}


template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const typename Aff_transformationC2<R CGAL_CTAG>::FT & m11, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m12,
                     const typename Aff_transformationC2<R CGAL_CTAG>::FT & m21, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m22,
                     const typename Aff_transformationC2<R CGAL_CTAG>::FT & w)
{
    PTR = new Aff_transformation_repC2<R>(m11/w, m12/w, m21/w, m22/w);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Translation, const typename Aff_transformationC2<R CGAL_CTAG>::Vector_2 &v)
{
  PTR = new Translation_repC2<R>(v);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2( const Rotation,
                      const typename Aff_transformationC2<R CGAL_CTAG>::Direction_2 &d,
                      const FT &num, const FT &den)
{
  PTR = new Rotation_repC2<R>(d, num, den);
}


template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Rotation,
                     const typename Aff_transformationC2<R CGAL_CTAG>::FT &sine, const typename Aff_transformationC2<R CGAL_CTAG>::FT &cosine, const typename Aff_transformationC2<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1)) // Idem...
    PTR = new Rotation_repC2<R>(sine/w, cosine/w);
  else
    PTR = new Rotation_repC2<R>(sine, cosine);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Scaling,
                     const typename Aff_transformationC2<R CGAL_CTAG>::FT &s, const typename Aff_transformationC2<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1)) // ....
    PTR = new Scaling_repC2<R>(s/w);
  else
    PTR = new Scaling_repC2<R>(s);
}



// and a 3x2 matrix for the operations combining rotation, scaling, translation
template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2( const typename Aff_transformationC2<R CGAL_CTAG>::FT & m11, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m12, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m13,
                      const typename Aff_transformationC2<R CGAL_CTAG>::FT & m21, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m22, const typename Aff_transformationC2<R CGAL_CTAG>::FT & m23,
                      const typename Aff_transformationC2<R CGAL_CTAG>::FT & w)
{
  if (w != FT(1)) // ...
    PTR = new Aff_transformation_repC2<R>(m11/w, m12/w, m13/w,
                                          m21/w, m22/w, m23/w);
  else
    PTR = new Aff_transformation_repC2<R>(m11, m12, m13,
                                          m21, m22, m23);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::~Aff_transformationC2()
{}

template < class R >
Aff_transformationC2<R CGAL_CTAG> &
Aff_transformationC2<R CGAL_CTAG>::operator=(const Aff_transformationC2<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
}

template < class R >
typename Aff_transformationC2<R CGAL_CTAG>::Point_2
Aff_transformationC2<R CGAL_CTAG>::transform(const typename Aff_transformationC2<R CGAL_CTAG>::Point_2 &p) const
{
  return ptr()->transform(p);
}

template < class R >
inline
typename Aff_transformationC2<R CGAL_CTAG>::Point_2
Aff_transformationC2<R CGAL_CTAG>::operator()(const typename Aff_transformationC2<R CGAL_CTAG>::Point_2 &p) const
{
  return transform(p);
}

template < class R >
typename Aff_transformationC2<R CGAL_CTAG>::Vector_2
Aff_transformationC2<R CGAL_CTAG>::transform(const typename Aff_transformationC2<R CGAL_CTAG>::Vector_2 &p) const
{
  return ptr()->transform(p);
}

template < class R >
inline
typename Aff_transformationC2<R CGAL_CTAG>::Vector_2
Aff_transformationC2<R CGAL_CTAG>::operator()(const typename Aff_transformationC2<R CGAL_CTAG>::Vector_2 &p) const
{
  return transform(p);
}
template < class R >
typename Aff_transformationC2<R CGAL_CTAG>::Direction_2
Aff_transformationC2<R CGAL_CTAG>::transform(const typename Aff_transformationC2<R CGAL_CTAG>::Direction_2 &d) const
{
  return ptr()->transform(d);
}

template < class R >
inline
typename Aff_transformationC2<R CGAL_CTAG>::Direction_2
Aff_transformationC2<R CGAL_CTAG>::operator()(const typename Aff_transformationC2<R CGAL_CTAG>::Direction_2 &d) const
{
  return transform(d);
}

template < class R >
inline
typename Aff_transformationC2<R CGAL_CTAG>::Line_2
Aff_transformationC2<R CGAL_CTAG>::transform(const typename Aff_transformationC2<R CGAL_CTAG>::Line_2 &l) const
{
  return typename Aff_transformationC2<R CGAL_CTAG>::Line_2(ptr()->transform(l.point(0)),
                    ptr()->transform(l.direction()));
}

template < class R >
inline
typename Aff_transformationC2<R CGAL_CTAG>::Line_2
Aff_transformationC2<R CGAL_CTAG>::operator()(const typename Aff_transformationC2<R CGAL_CTAG>::Line_2 &l) const
{
  return transform(l);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>
Aff_transformationC2<R CGAL_CTAG>::inverse() const
{
  return ptr()->inverse();
}


template < class R >
bool
Aff_transformationC2<R CGAL_CTAG>::is_odd() const
{
  return ! (ptr()->is_even());
}


template < class R >
std::ostream&
Aff_transformationC2<R CGAL_CTAG>::print(std::ostream &os) const
{
  ptr()->print(os);
  return os;
}


#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2
template < class R >
std::ostream&
operator<<(std::ostream& os, const Aff_transformationC2<R CGAL_CTAG>& t)
{
    t.print(os);
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONC2

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC2
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONC2

CGAL_END_NAMESPACE

#endif
