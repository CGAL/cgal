#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTCD_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTCD_H

#include <CGAL/number_utils.h>
#include <CGAL/determinant.h>
#include <CGAL/constructions/kernel_ftCd.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template < class InputIterator, class OutputIterator >
void
plane_from_pointsCd(int dim,
                    const InputIterator &db, const InputIterator &de,
                    OutputIterator result)
{
}

template < class InputIterator, class OutputIterator >
void
plane_from_point_directionCd(int dim,
                    const InputIterator &db1, const InputIterator &de1,
                    const InputIterator &db2, const InputIterator &de2,
                    OutputIterator result)
{
}

template < class InputIterator, class OutputIterator >
void
point_on_planeCd(int dim,
                 const InputIterator &db, const InputIterator &de,
                 int i, OutputIterator result)
{
}

template < class InputIterator, class OutputIterator >
void
projection_planeCd(int dim,
                   const InputIterator &hf, const InputIterator &hl,
                   const InputIterator &pf, const InputIterator &pl,
                   OutputIterator result)
{
}

CGAL_END_NAMESPACE

#endif  // CGAL_CONSTRUCTIONS_KERNEL_FTCD_H
