#ifndef CGAL_PREDICATES_KERNEL_FTCD_H
#define CGAL_PREDICATES_KERNEL_FTCD_H

#include <CGAL/number_utils.h>
#include <CGAL/determinant.h>
#include <CGAL/constructions/kernel_ftCd.h>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

template < class FT >
struct is_proportional { 
  FT _a, _b;
  is_proportional(const FT &a, const FT &b) : _a(a), _b(b) {}
  bool operator()(const FT &x1, const FT &x2) const
    { return (_a * x2 == _b * x1); }
};

template < class InputIterator >
CGAL_KERNEL_LARGE_INLINE
bool
is_positively_proportionalCd(
  const InputIterator &db1, const InputIterator &de1,
  const InputIterator &db2, const InputIterator &de2)
{
  typedef typename iterator_traits<InputIterator>::value_type FT;
  CGAL_kernel_assertion( de1-db1 == de2-db2 );
  InputIterator nz1; // first non-zero position in [db1,de1)
  nz1 = std::find_if(db1,de1,bind1st(not_equal_to<FT>(),FT(0)));
  CGAL_kernel_assertion( nz1 != de1 );
  InputIterator nz2 = db2 + (nz1-db1); // corresponding position in [db2,de2)
  // check that the factor of proportionality (*nz2/*nz1) is positive
  if ( *nz1 * *nz2 <= FT(0))
      return false;
  // check that all positions before nz2 are null 
  if (std::find_if(db2,nz2,bind1st(not_equal_to<FT>(),FT(0)))!=nz2)
      return false;
  // check that all subsequent positions are proportional
  return std::mismatch(nz1+1,de1+0,nz2+1,is_proportional<FT>(*nz1,*nz2))
	  .first == de1;
}

template < class InputIterator >
inline
bool
equal_directionCd(const InputIterator &db1, const InputIterator &de1,
                  const InputIterator &db2, const InputIterator &de2)
{
  return is_positively_proportionalCd(db1,de1,db2,de2);
}

CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_KERNEL_FTCD_H
