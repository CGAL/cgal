// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri, Lutz Kettner

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C

#include <CGAL/determinant.h>
#include <CGAL/Cartesian/Aff_transformation_rep_2.C>

#ifndef CGAL_CARTESIAN_REDEFINE_NAMES_2_H
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2()
{
  PTR = new Aff_transformation_repC2<R>(FT(1), FT(0), FT(0), FT(1));
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Identity_transformation)
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
Aff_transformationC2(
        const typename Aff_transformationC2<R CGAL_CTAG>::FT & m11,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT & m12,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT & m21,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT & m22,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT & w)
{
    PTR = new Aff_transformation_repC2<R>(m11/w, m12/w, m21/w, m22/w);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(
	const Translation,
        const typename Aff_transformationC2<R CGAL_CTAG>::Vector_2 &v)
{
  PTR = new Translation_repC2<R>(v);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(
        const Rotation,
        const typename Aff_transformationC2<R CGAL_CTAG>::Direction_2 &d,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT &num,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &den)
{
  PTR = new Rotation_repC2<R>(d, num, den);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(
        const Rotation,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT &sine,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &cosine,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1)) // Idem...
    PTR = new Rotation_repC2<R>(sine/w, cosine/w);
  else
    PTR = new Rotation_repC2<R>(sine, cosine);
}

template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(const Scaling,
                     const typename Aff_transformationC2<R CGAL_CTAG>::FT &s,
		     const typename Aff_transformationC2<R CGAL_CTAG>::FT &w)
{
  if (w != FT(1)) // ....
    PTR = new Scaling_repC2<R>(s/w);
  else
    PTR = new Scaling_repC2<R>(s);
}

// and a 3x2 matrix for the operations combining rotation, scaling, translation
template < class R >
Aff_transformationC2<R CGAL_CTAG>::
Aff_transformationC2(
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &m11,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT &m12,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &m13,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT &m21,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &m22,
	const typename Aff_transformationC2<R CGAL_CTAG>::FT &m23,
        const typename Aff_transformationC2<R CGAL_CTAG>::FT &w)
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
Aff_transformationC2<R CGAL_CTAG>::
operator=(const Aff_transformationC2<R CGAL_CTAG> &t)
{
  Handle::operator=(t);
  return *this;
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

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_2_C
