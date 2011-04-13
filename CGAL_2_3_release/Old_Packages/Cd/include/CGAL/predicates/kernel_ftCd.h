// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/predicates/kernel_ftCd.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_PREDICATES_KERNEL_FTCD_H
#define CGAL_PREDICATES_KERNEL_FTCD_H

#include <CGAL/number_utils.h>
#include <CGAL/determinant.h>
#include <CGAL/constructions/kernel_ftCd.h>
#include <algorithm>
#include <numeric>

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
is_proportionalCd(const InputIterator &first1, const InputIterator &last1,
                  const InputIterator &first2)
{
  typedef typename iterator_traits<InputIterator>::value_type FT;
  InputIterator nz1; // first non-zero position in [first1,last1)
  nz1 = std::find_if(first1,last1,bind1st(not_equal_to<FT>(),FT(0)));
  CGAL_kernel_assertion( nz1 != last1 );
  InputIterator nz2 = first2 + (nz1-first1);
  // check that all positions before nz2 are null 
  if (std::find_if(first2,nz2,bind1st(not_equal_to<FT>(),FT(0)))==nz2)
  // check that all subsequent positions are proportional
    return std::mismatch(nz1+1,last1+0,nz2+1,
                         is_proportional<FT>(*nz1,*nz2)).first == last1;
  // otherwise
  return false;
}

template < class InputIterator >
CGAL_KERNEL_LARGE_INLINE
bool
is_positively_proportionalCd(
  const InputIterator &first1, const InputIterator &last1,
  const InputIterator &first2)
{
  typedef typename iterator_traits<InputIterator>::value_type FT;
  InputIterator nz1; // first non-zero position in [first1,last1)
  nz1 = std::find_if(first1,last1,bind1st(not_equal_to<FT>(),FT(0)));
  CGAL_kernel_assertion( nz1 != last1 );
  InputIterator nz2 = first2 + (nz1-first1);
  // check that the factor of proportionality (*nz2/*nz1) is positive
  if ( CGAL_NTS sign(*nz1) == CGAL_NTS sign(*nz2))
  // check that all positions before nz2 are null 
    if (std::find_if(first2,nz2,bind1st(not_equal_to<FT>(),FT(0)))==nz2)
  // check that all subsequent positions are proportional
      return std::mismatch(nz1+1,last1+0,nz2+1,
                           is_proportional<FT>(*nz1,*nz2)).first == last1;
  // otherwise
  return false;
}

template < class InputIterator >
Comparison_result
compare_dominanceCd(const InputIterator &/*pb*/, const InputIterator &/*pe*/,
                    const InputIterator &/*qb*/, const InputIterator &/*qe*/)
{
  // TODO
  return SMALLER;
}

template < class InputIterator >
Comparison_result
compare_submittanceCd(const InputIterator &/*pb*/, const InputIterator &/*pe*/,
                      const InputIterator &/*qb*/, const InputIterator &/*qe*/)
{
  // TODO
  return SMALLER;
}

template < class R, class PointIterator >
Orientation
orientationCd(const PointIterator &first, const PointIterator &last,
              const R&)
{
  typename R::LA::Matrix M(first->dimension());
  // subtract first point to each point and put it in consecutive rows
  typename R::LA::Matrix::const_iterator row;
  PointIterator p; 
  for (row=M.begin(), p=first+1; p!=last; ++p) {
    CGAL_kernel_precondition( p->dimension() == first->dimension() );
    row = std::transform(p->begin(), p->end(), first->begin(),
                         row, std::minus<typename R::FT>());
  }
  // orientation is the sign of the determinant
  typename R::LA La;
  return Orientation(La.sign_of_determinant(M));
}

template < class InputIterator >
bool
collinear_are_ordered_along_lineCd(
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &rb)
{
  InputIterator p, q, r;
  for (p=pb,q=qb,r=rb; p!=pe; ++p,++q,++r) {
    if (*p<*q) return !(*r<*q);
    if (*q<*p) return !(*q<*r);
  }
  return true; // p==q
}

template < class InputIterator >
bool
collinear_are_strictly_ordered_along_lineCd(
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &rb)
{
  InputIterator p, q, r;
  for (p=pb,q=qb,r=rb; p!=pe; ++p,++q,++r) {
    if (*p<*q) return *q<*r;
    if (*q<*p) return *r<*q;
  }
  return false; // p==q
}

template < class InputIterator >
inline
bool
equal_directionCd(const InputIterator &db1, const InputIterator &de1,
                  const InputIterator &db2)
{
  return is_positively_proportionalCd(db1,de1,db2);
}

template < class InputIterator >
Oriented_side
side_of_oriented_planeCd(const InputIterator &hb, const InputIterator &he,
                         const InputIterator &pb)
{
  std::iterator_traits<InputIterator>::value_type zero(0);
  return Oriented_side(
            CGAL_NTS compare(std::inner_product(hb, he-1, pb, *(he-1)), zero) );
}

template < class R, class PointIterator, class InputIterator >
Bounded_side
side_of_bounded_simplexCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,
   const R&)
{
  typedef typename R::FT FT;
  typename R::LA::Matrix M(first->dimension());
  typename R::LA::Vector b(first->dimension());
  // subtract last point to each point and put it in consecutive rows
  typename R::LA::Matrix::const_iterator row;
  PointIterator p, lastp = last-1; 
  for (row=M.begin(), p=first; p!=lastp; ++p) {
    CGAL_kernel_precondition( p->dimension() == first->dimension() );
    row = std::transform(p->begin(), p->end(), lastp->begin(),
                         row, std::minus<FT>());
  }
  // substract last point to query point and put it in b
  CGAL_kernel_assertion( pe-pb == first->dimension() );
  std::transform(pb, pe, lastp->begin(), b.begin(), std::minus<FT>());
  // affine coordinates are solutions of the system Mx=b
  typename R::LA::Vector x(first->dimension());
  typename R::LA La;
  FT D;
  bool is_solvable = La.linear_solver(La.transpose(M),b,x,D);
  CGAL_kernel_assertion( is_solvable );
  if (D<FT(0)) {
    D = -D;
    std::transform(x.begin(),x.end(),x.begin(),std::negate<FT>());
  }
  bool degenerate = false;
  typename R::LA::Vector::const_iterator xit;
  for (xit=x.begin(); xit!=x.end(); ++xit) {
    if (*xit < FT(0)) return ON_UNBOUNDED_SIDE;
    if (*xit == FT(0)) degenerate = true;
  }
  if (degenerate) return ON_BOUNDARY;
  FT sum = std::accumulate(x.begin(),x.end(),FT(0));
  if (sum > D) return ON_UNBOUNDED_SIDE;
  if (sum == D) return ON_BOUNDARY;
  return ON_BOUNDED_SIDE;
}

template < class R, class PointIterator, class InputIterator >
Oriented_side
side_of_oriented_sphereCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,
   const R&)
{
  typedef typename R::FT FT;
  typename R::LA::Matrix M(first->dimension()+1);
  // subtract last point to each point and put it in consecutive rows
  // put sum of squares at the end of each column
  PointIterator p, lastp = last-1;
  typename R::LA::Matrix::const_iterator row, row2;
  for (row=M.begin(), p=first; p!=lastp; ++p) {
    row2 = std::transform(p->begin(),p->end(),lastp->begin(),
                          row,std::minus<typename R::FT>()); //? ++row2;
    *row2 = std::inner_product(row,row+first->dimension(),row,FT(0)); ++row2;
    row = row2;
  }
  // substract last point to query point
  // and put sum of squares at the end of each column
  CGAL_kernel_assertion( pe-pb == first->dimension() );
  row2 = std::transform(pb, pe, lastp->begin(),row,std::minus<FT>());
  *row2 = std::inner_product(row,row+first->dimension(),row,FT(0));
  // side_of_sphere is the sign of the determinant of M
  typename R::LA La;
  return Oriented_side(La.sign_of_determinant(M));
}

template < class R, class PointIterator, class InputIterator >
Bounded_side
side_of_bounded_sphereCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,   const R&r)
{
  Oriented_side s = side_of_oriented_sphereCd(first,last,pb,pe,r);
  Orientation o = orientationCd(first,last,r);
  return Bounded_side(s * o);
}

template < class InputIterator >
Comparison_result
cmp_dist_to_pointCd(
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &/*qe*/,
   const InputIterator &rb, const InputIterator &/*re*/)
{
  return CGAL_NTS compare(squared_distanceCd(pb,pe,qb),
                          squared_distanceCd(pb,pe,rb));
}

template < class InputIterator >
bool
has_larger_dist_to_pointCd(
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe,
   const InputIterator &rb, const InputIterator &re)
{
  return cmp_dist_to_pointCd(pb,pe,qb,qe,rb,re) == LARGER;
}

template < class InputIterator >
bool
has_smaller_dist_to_pointCd(
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe,
   const InputIterator &rb, const InputIterator &re)
{
  return cmp_dist_to_pointCd(pb,pe,qb,qe,rb,re) == SMALLER;
}

template < class InputIterator >
Comparison_result
cmp_signed_dist_to_directionCd(
   const InputIterator &hb, const InputIterator &he,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return CGAL_NTS compare(scaled_distance_to_directionCd(hb,he,pb,pe),
                          scaled_distance_to_directionCd(hb,he,qb,qe));
}

template < class InputIterator >
bool
has_larger_signed_dist_to_directionCd(
   const InputIterator &hb, const InputIterator &he,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return cmp_signed_dist_to_directionCd(hb,he,pb,pe,qb,qe) == LARGER;
}

template < class InputIterator >
bool
has_smaller_signed_dist_to_directionCd(
   const InputIterator &hb, const InputIterator &he,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return cmp_signed_dist_to_directionCd(hb,he,pb,pe,qb,qe) == SMALLER;
}

template < class PointIterator, class InputIterator >
Comparison_result
cmp_signed_dist_to_directionCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return CGAL_NTS compare(scaled_distance_to_directionCd(first,last,pb,pe),
                          scaled_distance_to_directionCd(first,last,qb,qe));
}

template < class PointIterator, class InputIterator >
bool
has_larger_signed_dist_to_directionCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return cmp_signed_dist_to_directionCd(first,last,pb,pe,qb,qe) == LARGER;
}

template < class PointIterator, class InputIterator >
bool
has_smaller_signed_dist_to_directionCd(
   const PointIterator &first, const PointIterator &last,
   const InputIterator &pb, const InputIterator &pe,
   const InputIterator &qb, const InputIterator &qe)
{
  return cmp_signed_dist_to_directionCd(first,last,pb,pe,qb,qe) == SMALLER;
}

CGAL_END_NAMESPACE

#endif  // CGAL_PREDICATES_KERNEL_FTCD_H
