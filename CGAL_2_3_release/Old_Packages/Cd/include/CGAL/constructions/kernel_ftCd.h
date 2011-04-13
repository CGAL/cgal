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
// file          : include/CGAL/constructions/kernel_ftCd.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Hervé Brönnimann
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

#ifndef CGAL_CONSTRUCTIONS_KERNEL_FTCD_H
#define CGAL_CONSTRUCTIONS_KERNEL_FTCD_H

#include <CGAL/number_utils.h>
#include <CGAL/determinant.h>
#include <CGAL/constructions/kernel_ftCd.h>
#include <algorithm>
#include <functional>

CGAL_BEGIN_NAMESPACE

template < class R, class InputIterator, class OutputIterator >
OutputIterator
point_on_lineCd(const InputIterator &pb, const InputIterator &pe,
                const InputIterator &db, typename R::FT ratio,
	        OutputIterator result, const R&)
{
  // METHOD: result is simply lp + i*ld
  InputIterator ip, id;
  for (ip=pb,id=db; ip!=pe; ++ip,++id) {
    *result = *ip + ratio * (*id); ++result;
  }
  return result;
}

template < class R, class InputIterator, class OutputIterator >
OutputIterator
projection_lineCd(const InputIterator &pb,  const InputIterator &pe,
                  const InputIterator &lpb, const InputIterator &lpe,
                  const InputIterator &ldb, const InputIterator &lde,
	          OutputIterator result, const R &r)
{

  // METHOD: scalar products ratio ((p-lp)*ld/ld*ld) gives the
  // projection on the line in ld units, starting from lp
  // hence x = lp + ((p-lp)*ld/ld*ld) * ld
  typedef typename R::FT FT;
  InputIterator ip, ilp, ild;
  FT ratio( FT(0) );
  for (ip=pb,ilp=lpb,ild=ldb; ip!=pe; ++ip,++ilp,++ild)
    ratio += (*ip - *ilp) * (*ild);
  ratio /= std::inner_product(ldb,lde,ldb,FT(0));
  return point_on_lineCd(lpb, lpe, ldb, ratio, result, r);
}

template < class R, class PointIterator, class OutputIterator >
OutputIterator
plane_from_pointsCd(int dim,
   const PointIterator &pb, const PointIterator &pe, // points
   OutputIterator result, const R&)
{
  // METHOD: Write equation of plane as det(M')=0, where M is the matrix
  // obtained by writing the coordinates of the points in [pb..pe) on
  // each row, and M' obtained by adding on the last row the
  // indeterminates, and all 1s in the last column. The coefficient of
  // x[i], its minor, is also the determinant of M where the ith column
  // is replaced by all 1s; note that the (-1)^(d+i) is cancelled by the
  // (d-i) permutations needed to bring the 1s in the minor into the ith
  // column. The constant term corresponding to i==d, is the
  // determinant of M itself.
  CGAL_kernel_precondition( pe-pb == dim );
  typedef typename R::FT FT;
  // InputIterator::value_type is Point_d
  // OutputIterator::value_type is FT
  typename R::LA La;
  typename R::LA::Vector save(dim);
  typename R::LA::Matrix M(dim, pb, pe);
  int j;
  for (j=0; j<dim; ++j) {
    // Replace column(j) by (1...1), and save previous column
    int k;
    for (k=0; k!=dim; ++k) {
       save[k] = M[k][j];
       M[k][j] = FT(1);
    }
    // Compute minor corresponding to variable x_i and save into result
    *result = La.determinant(M); ++result;
    // Restore column
    for (k=0; k!=dim; ++k)
       M[k][j] = save[k];
  }
  // Compute determinant and save into result
  *result = La.determinant(M); ++result;
  return result;
}

template < class R, class InputIterator, class OutputIterator >
OutputIterator
plane_from_point_directionCd(int dim,
   const InputIterator &pb, const InputIterator &pe, // point
   const InputIterator &db, const InputIterator &de, // direction
   OutputIterator result, const R&)
{
  // METHOD: the first d coefficients are simply those of the direction,
  // the remaining one is such that the equation cancels for p.
  CGAL_kernel_precondition( pe-pb == dim );
  CGAL_kernel_precondition( de-db == dim );
  result = std::copy(db, de, result);
  *result = - std::inner_product(pb, pe, db, typename R::FT(0)); ++result;
  return result;
}

template < class R, class InputIterator, class OutputIterator >
OutputIterator
point_on_planeCd(int dim,
   const InputIterator &hb, const InputIterator &he, int i, 
   OutputIterator result, const R&)
{
  CGAL_kernel_precondition( he-hb == dim+1 );
  i = (i%dim); if (i<0) i += dim;
  // METHOD: if h[i]!=0, then the result is x[i]=-h[d]/h[i], x[j]=0
  // and if h[i]==0, we look for the first k such that h[k]!=0, and put
  // x[k]=-h[d]/h[k], x[i]=1, and x[j]=0 otherwise.
  typedef typename std::iterator_traits<InputIterator>::value_type FT;
  int j;
  if (*(hb+i)!=FT(0)) {
    for (j=0; j<i; ++j) { *result = FT(0); ++result; }
    *result = - *(hb+dim) / *(hb+i); ++result;
    for (j=i+1; j<dim; ++j) { *result = FT(0); ++result; }
  } else {
    InputIterator hk = std::find_if(hb,he-1,bind1st(not_equal_to<FT>(),FT(0)));
    CGAL_kernel_assertion( hk!=he-1 );
    if (i<hk-hb) {
      for (j=0; j<i; ++j) { *result = FT(0); ++result; }
      *result = FT(1); ++result;
      for (++j; j<hk-hb; ++j) { *result = FT(0); ++result; }
      *result = - *(hb+dim) / *hk; ++result;
      for (++j; j<dim; ++j) { *result = FT(0); ++result; }
    } else {
      for (j=0; j<hk-hb; ++j) { *result = FT(0); ++result; }
      *result = - *(hb+dim) / *hk; ++result;
      for (++j; j<i; ++j) { *result = FT(0); ++result; }
      *result = FT(1); ++result;
      for (++j; j<dim; ++j) { *result = FT(0); ++result; }
    }
  }
  return result;
}

template < class R, class InputIterator, class OutputIterator >
OutputIterator
projection_planeCd(int dim,
   const InputIterator &hb, const InputIterator &he, // plane
   const InputIterator &pb, const InputIterator &pe, // point
   OutputIterator result, const R&)
{
  // Look for x = x[0]...x[d-1] such that for i ranging in [0,d)
  // (1) sum_i(h[i]*x[i]) + h[d] = 0, and (2) x-p is collinear to h[0]...h[d-1]
  // METHOD: pick a coordinate for which h[k]!=FT(0),then
  // collinearity         (2) =>    x[i]-p[i] = (x[k]-p[k])*h[i]/h[k].
  // hyperplane eqn for p (3) =>    sum_i(h[i]*p[i]) + hp = 0, for some hp 
  // hyperplane eqn (1) - (3) =>    x[k]-p[k] = h[k]*(hp-h[d])/sum_k(h[i]^2)
  // hence x[i] = p[i] + h[i] * (hp-h[d]) / sum_k(h[i]^2), for any i
  typedef typename std::iterator_traits<InputIterator>::value_type FT;
  FT sum_of_sqr_h = std::inner_product(hb,he-1,hb,FT(0));
  FT hpmhd = std::inner_product(pb,pe,hb, *(hb+dim));
  int i;
  for (i=0; i<dim; ++i) {
    *result = *(pb+i) - *(hb+i) * hpmhd / sum_of_sqr_h; ++result;
  }
  return result;
}

template < class InputIterator >
typename std::iterator_traits<InputIterator>::value_type
squared_distanceCd(InputIterator pb, const InputIterator &pe,
                   InputIterator qb)
{
  typename std::iterator_traits<InputIterator>::value_type result(0);
  for (; pb!=pe; ++pb,++qb)
    result += CGAL_NTS square(*pb - *qb);
  return result;
}

template < class InputIterator >
typename std::iterator_traits<InputIterator>::value_type
scaled_distance_to_directionCd(const InputIterator &hb,
                               const InputIterator &pb,
			       const InputIterator &pe)
{
  typedef typename std::iterator_traits<InputIterator>::value_type FT;
  return std::inner_product (pb,pe,hb,FT(0));
}

template < class InputIterator >
typename std::iterator_traits<InputIterator>::value_type
scaled_distance_to_planeCd(const InputIterator &hb, const InputIterator &he,
                           const InputIterator &pb, const InputIterator &pe)
{
  return std::inner_product (pb,pe,hb,*(he-1));
}

template < class R, class PointIterator, class InputIterator >
typename std::iterator_traits<InputIterator>::value_type
scaled_distance_to_planeCd(const PointIterator &first,
                           const PointIterator &last,
                           const InputIterator &pb, const InputIterator &pe,
			   const R&)
{
  typename R::LA::Matrix M(first->dimension(), first, last);
  // subtract point [pb,pe) to each row
  typename R::LA::Matrix::const_iterator row;
  for (row=M.begin(); row!=M.end(); )
    row = std::transform(row, row+first->dimension(), pb,
                         row, std::minus<typename R::FT>());
  // scaled distance is simply the determinant
  typename R::LA La;
  return La.determinant(M);
}

CGAL_END_NAMESPACE

#endif  // CGAL_CONSTRUCTIONS_KERNEL_FTCD_H
