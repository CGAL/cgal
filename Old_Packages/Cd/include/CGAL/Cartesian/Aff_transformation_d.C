// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Brönnimann

#ifndef CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C
#define CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C

#include <CGAL/Cartesian/Aff_transformation_rep_d.C>
#include <CGAL/determinant.h>

#ifndef CGAL_CTAG
#define CGAL_CTAG
#endif

#ifdef CGAL_CFG_TYPENAME_BUG
#define typename
#endif

CGAL_BEGIN_NAMESPACE

template < class R >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd()
{
}

template < class R >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd(const Identity_transformation &)
{
}

template < class R >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd(
	const Translation,
        const typename Aff_transformationCd<R CGAL_CTAG>::Vector_d &v)
{
}

template < class R >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd(const Scaling,
                     const typename Aff_transformationCd<R CGAL_CTAG>::FT &s,
                     const typename Aff_transformationCd<R CGAL_CTAG>::FT &w)
{
}

template < class R, class InputIterator >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd(int d,
                     const InputIterator &begin, const InputIterator &last,
                     const typename Aff_transformationCd<R CGAL_CTAG>::FT &w)
{
}

template < class R, class InputIterator >
Aff_transformationCd<R CGAL_CTAG>::
Aff_transformationCd(int d,
                     const InputIterator &begin, const InputIterator &last,
                     const InputIterator &translation_begin,
		     const InputIterator &translation_last,
                     const typename Aff_transformationCd<R CGAL_CTAG>::FT &w)
{
}

template < class R >
Aff_transformationCd<R CGAL_CTAG>::~Aff_transformationCd()
{}

#ifndef CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONCD
template < class R >
std::ostream &operator<<(std::ostream &os,
                         const Aff_transformationCd<R CGAL_CTAG> &t)
{
    t.print(os);
    return os;
}
#endif // CGAL_NO_OSTREAM_INSERT_AFF_TRANSFORMATIONCD

#ifndef CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONCD
#endif // CGAL_NO_ISTREAM_EXTRACT_AFF_TRANSFORMATIONCD

CGAL_END_NAMESPACE

#ifdef CGAL_CFG_TYPENAME_BUG
#undef typename
#endif

#endif // CGAL_CARTESIAN_AFF_TRANSFORMATION_D_C
