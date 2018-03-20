// Copyright (c) 2000,2001  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
// 
//
// Author(s)     : Michael Seel

#ifndef CGAL_PREDICATES_D_H
#define CGAL_PREDICATES_D_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>
#include <CGAL/Kernel_d/Point_d.h>

namespace CGAL {

/*{\Moptions outfile=predicates_d.man}*/
/*{\Mtext \setopdims{4cm}{2cm}\computewidths
\headerline{Linear and affine predicates}

For an iterator range |[first,last)| we define |T = tuple [first,last)|
as the ordered tuple $(|T[0]|,|T[1]|, \ldots |T[d-1]|)$ where 
$|T[i]| = |*| |++|^{(i)}|first|$ (the element obtained by $i$ times 
forwarding the iterator by operator |++| and then dereferencing it to
get the value to which it points). We write |d = size [first,last)|.
This extends the syntax of random access iterators to input iterators. 
If we index the tuple as above then we require that 
$|++|^{(d+1)}|first == last|$.

In the following we require the Iterators to be globally 
of value type |Point_d<R>|. Also if we are handed over an iterator
range |[first,last)|, then all points in |T = tuple [first,last)|
have the same dimension |dim|.
}*/

template <class R, class ForwardIterator, class OutputIterator>
OutputIterator barycentric_coordinates(
  ForwardIterator first, ForwardIterator last, const Point_d<R>& p, 
  OutputIterator result)
/*{\Xfunc returns the barycentric coordinates of a point $p \in R^d$ in a
affine space of dimension $k$ spanned by the points in |tuple [first,last)|.
\precond value type of |ForwardIterator| is |Point_d<R>|,
|affinely_independed(first,last)| and 
|affine_rank(tuple [first,last),p)==k|.}*/
{ typename R::Barycentric_coordinates_d coords;
  return coords(first,last,p,result);
}

template <class ForwardIterator>
Orientation
orientation(ForwardIterator first, ForwardIterator last)
/*{\Mfunc determines the orientation of the points in the tuple 
|A = tuple [first,last)| where $A$ consists of $d + 1$ points in $d$ - space. 
This is the sign of the determinant
  \[ \left\Lvert \begin{array}{cccc}
  1 & 1 & 1 & 1 \\
  A[0] & A[1] & \dots & A[d]
  \end{array}  \right\Lvert  \]
  where |A[i]| denotes the cartesian coordinate vector of 
  the $i$-th point in $A$.
  \precond |size [first,last) == d+1| and |A[i].dimension() == d| 
  $\forall 0 \leq i \leq d$
  and the value type of |ForwardIterator| is |Point_d<R>|.}*/
{ typedef typename std::iterator_traits<ForwardIterator>::value_type PointD;
  typedef typename PointD::R R;
  typename R::Orientation_d or_; return or_(first,last); }


template <class ForwardIterator>
Orientation
coaffine_orientation(ForwardIterator first, ForwardIterator last)
/*{\Mfunc determines the orientation of the points in the tuple 
|A = tuple [first,last)| where $A$ consists of $k + 1$ points in $d$ - space. 
This is the sign of the determinant
  \[ \left\Lvert \begin{array}{cccc}
  1 & 1 & 1 & 1 \\
  A[0] & A[1] & \dots & A[k]
  \end{array}  \right\Lvert  \]
  where |A[a_i]| denotes the cartesian coordinate vector of 
  the $i$-th point in $A$, after projection on an axis-aligned |k|-affine subspace.
  \precond |size [first,last) <= d| and |A[i].dimension() == d| 
  $\forall 0 \leq i \leq d$
  and the value type of |ForwardIterator| is |Point_d<R>|.}*/
{
    typedef typename std::iterator_traits<ForwardIterator>::value_type PointD;
    typedef typename PointD::R R;
    typename R::Coaffine_orientation_d or_;
    return or_(first,last);
}


template <class R, class ForwardIterator>
Oriented_side side_of_oriented_sphere(
  ForwardIterator first, ForwardIterator last, 
  const Point_d<R>& x)
/*{\Mfuncl determines to which side of the sphere $S$ defined by the
    points in |A = tuple [first,last)| the point $x$ belongs, 
    where $A$ consists of $d + 1$ points in $d$ - space.  
    The positive side is determined by the positive sign 
    of the determinant \[
    \left\Lvert \begin{array}{ccccc} 
    1 & 1 & 1 & 1 & 1\\
    |lift(A[0])| & |lift(A[1])| & \dots & |lift(A[d])| & |lift(x)| 
    \end{array} \right\Lvert \] 
    where for a point $p$ with cartesian coordinates $p_i$ we use 
    |lift(p)| to denote the $d + 1$-dimensional point with cartesian 
    coordinate vector $(p_0, \ldots,p_{d-1},\sum_{0 \le i < d}p_i^2)$. If
    the points in $A$ are positively oriented then the positive
    side is the inside of the sphere and the negative side is
    the outside of the sphere.  
\precond value type of |ForwardIterator| is |Point_d<R>|.}*/
{ typename R::Side_of_oriented_sphere_d side; 
  return side(first,last,x); }

template <class R, class ForwardIterator>
Bounded_side side_of_bounded_sphere(
  ForwardIterator first, ForwardIterator last,
  const Point_d<R> &p)
/*{\Mfunc determines whether the point |p| lies |ON_BOUNDED_SIDE|, 
|ON_BOUNDARY|, or |ON_UNBOUNDED_SIDE| of the sphere defined by 
the points in |A = tuple [first,last)| where $A$ consists of $d+1$ 
points in $d$-space.
\precond value type of |ForwardIterator| is |Point_d<R>| and
$|orientation(first,last)| \neq |ZERO|$.}*/
{ typename R::Side_of_bounded_sphere_d side; 
  return side(first,last,p); }

template <class R, class ForwardIterator>
bool contained_in_simplex(
  ForwardIterator first, ForwardIterator last, const Point_d<R>& p) 
/*{\Mfunc determines whether |p| is contained in the simplex spanned
by the points in |A = tuple [first,last)|. |A| may consists of up to 
$d + 1$ points. 
\precond value type of |ForwardIterator| is |Point_d<R>| and
the points in |A| are affinely independent.}*/
{ typename R::Contained_in_simplex_d contained; 
  return contained(first,last,p); }

template <class R, class ForwardIterator>
bool contained_in_affine_hull(
     ForwardIterator first, ForwardIterator last,
     const Point_d<R>& p) 
/*{\Mfunc determines whether $p$ is contained in the affine hull
of the points in |A = tuple [first,last)|.
\precond value type of |ForwardIterator| is |Point_d<R>|.}*/
{ typename R::Contained_in_affine_hull_d contained; 
  return contained(first,last,p); }

template <class ForwardIterator>
int affine_rank(ForwardIterator first, ForwardIterator last)
/*{\Mfunc computes the affine rank of the points in |A = tuple [first,last)|.
\precond value type of |ForwardIterator| is |Point_d<R>|.}*/
{ typedef typename std::iterator_traits<ForwardIterator>::value_type PointD;
  typedef typename PointD::R R;
  typename R::Affine_rank_d rank; 
  return rank(first,last); }

template <class ForwardIterator>
bool affinely_independent(ForwardIterator first, ForwardIterator last)
/*{\Mfunc decides whether the points in |A = tuple [first,last)| are 
affinely independent. 
\precond value type of |ForwardIterator| is |Point_d<R>|.}*/
{ typedef typename std::iterator_traits<ForwardIterator>::value_type PointD;
  typedef typename PointD::R R;
  typename R::Affinely_independent_d ind; 
  return ind(first,last); }

template <class R>
Comparison_result compare_lexicographically( 
  const Point_d<R>& p1, const Point_d<R>& p2)
/*{\Mfunc compares the Cartesian coordiantes of points |p1| and |p2|
   lexicographically.}*/
{ typename R::Compare_lexicographically_d cmp; 
  return cmp(p1,p2); }

template <class R>
bool lexicographically_smaller( 
  const Point_d<R>& p1, const Point_d<R>& p2)
/*{\Mfunc returns true iff $|p1| < |p2|$ in the Cartesian lexicographic
   order of points.}*/
{ typename R::Less_lexicographically_d lt;
  return lt(p1,p2); }

template <class R>
bool lexicographically_smaller_or_equal(
  const Point_d<R>& p1, const Point_d<R>& p2)
/*{\Mfunc returns true iff $|p1| \leq |p2|$ in the Cartesian lexicographic
   order of points.}*/
{ typename R::Less_or_equal_lexicographically_d le;
  return le(p1,p2); }

template <class R, class ForwardIterator>
bool contained_in_linear_hull(
  ForwardIterator first, ForwardIterator last, const Vector_d<R>& x)
/*{\Mfunc determines whether $x$ is contained in the linear hull
of the vectors in |A = tuple [first,last)|. 
\precond value type of |ForwardIterator| is |Vector_d<R>|.}*/
{ typename R::Contained_in_linear_hull_d contained; 
  return contained(first,last,x); }

template <class ForwardIterator>
int linear_rank(
  ForwardIterator first, ForwardIterator last)
/*{\Mfunc computes the linear rank of the vectors in 
|A = tuple [first,last)|. 
\precond value type of |ForwardIterator| is |Vector_d<R>|.}*/
{ typedef typename std::iterator_traits<ForwardIterator>::
    value_type value_type;
  typedef typename value_type::R R;
  typename R::Linear_rank_d rank; 
  return rank(first,last); 
}

template <class ForwardIterator>
bool linearly_independent(
  ForwardIterator first, ForwardIterator last)
/*{\Mfunc decides whether the vectors in $A$ are linearly independent.
\precond value type of |ForwardIterator| is |Vector_d<R>|.}*/
{ TUPLE_DIM_CHECK(first,last,linearly_independent);
  typedef typename std::iterator_traits<ForwardIterator>::
    value_type value_type;
  typedef typename value_type::R R;
  typename R::Linearly_independent_d independent; 
  return independent(first,last);
}


} //namespace CGAL
#endif // CGAL_PREDICATES_D_H
