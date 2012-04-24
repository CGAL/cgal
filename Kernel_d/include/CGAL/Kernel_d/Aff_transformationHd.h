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
// 
//
// Author(s)     : Michael Seel

#ifndef CGAL_AFF_TRANSFORMATIONHD_H
#define CGAL_AFF_TRANSFORMATIONHD_H

#include <CGAL/basic.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/rational_rotation.h>
#include <CGAL/Handle_for.h>

namespace CGAL {

template <class RT, class LA > class Aff_transformationHd;
template <class RT, class LA > class Aff_transformationHd_rep;

template <class RT, class LA>
class Aff_transformationHd_rep 
{
  friend class Aff_transformationHd<RT,LA>;
  typedef typename LA::Matrix Matrix;
  Matrix M_;
public:
  Aff_transformationHd_rep(int d) : M_(d+1) {}
  Aff_transformationHd_rep(const Matrix& M_init) : M_(M_init) {}
  ~Aff_transformationHd_rep() {}
}; 


/*{\Moptions outfile=Aff_transformation_d.man}*/ 
/*{\Manpage{Aff_transformation_d}{R}{Affine Transformations}{t}}*/
/*{\Msubst 
Hd<RT,LA>#_d<R>
Aff_transformationHd#Aff_transformation_d
Quotient<RT>#FT
}*/

template <class _RT, class _LA>
class Aff_transformationHd : 
  public Handle_for< Aff_transformationHd_rep<_RT,_LA> > { 

  typedef Aff_transformationHd_rep<_RT,_LA> Rep;
  typedef Handle_for<Rep> Base;
  typedef Aff_transformationHd<_RT,_LA> Self;

  using Base::ptr;

/*{\Mdefinition 
An instance of the data type |\Mname| is an affine transformation of
$d$-dimensional space. It is specified by a square matrix
$M$ of dimension $d + 1$. All entries in the last row of |M| except
the diagonal entry must be zero; the diagonal entry must be non-zero.
A point $p$ with homogeneous coordinates $(p[0], \ldots, p[d])$ can be
transformed into the point |p.transform(A)|, where |A| is an affine
transformation created from |M| by the constructors below. }*/

public: 
/*{\Mtypes 4}*/

typedef _RT RT;
/*{\Mtypemember the ring type.}*/
typedef Quotient<_RT> FT;
/*{\Mtypemember the field type.}*/
typedef _LA LA;
/*{\Mtypemember the linear algebra layer.}*/
typedef typename _LA::Matrix Matrix;
/*{\Mtypemember the matrix type.}*/
typedef typename _LA::Vector Vector;

/*{\Mcreation 3}*/

Aff_transformationHd(int d = 0) : Base( Rep(d) ) {}
/*{\Mcreate introduces a transformation in $d$-dimensional space.}*/

Aff_transformationHd(int d, Identity_transformation) : Base( Rep(d) )
/*{\Mcreate introduces the identity transformation in $d$-dimensional 
    space.}*/
{ for (int i = 0; i <= d; ++i) ptr()->M_(i,i) = RT(1); }

Aff_transformationHd(const Matrix& M) : Base( Rep(M) )
/*{\Mcreate introduces the transformation of $d$ - space specified by
matrix $M$. \precond |M| is a square matrix of dimension $d + 1$. }*/
{ CGAL_assertion_msg((M.row_dimension()==M.column_dimension()),
    "Aff_transformationHd::\
     construction: initialization matrix is not quadratic.");
}

template <typename Forward_iterator>
Aff_transformationHd(Scaling, Forward_iterator start, Forward_iterator end) :
  Base( Rep(static_cast<int>(std::distance(start,end))-1) )
/*{\Mcreate introduces the transformation of $d$-space specified by a
diagonal matrix with entries |set [start,end)| on the diagonal 
(a scaling of the space). \precond |set [start,end)| is a vector of 
dimension $d+1$.}*/
{ int i=0; while (start != end) { ptr()->M_(i,i) = *start++;++i; } }

Aff_transformationHd(Translation, const VectorHd<RT,LA>& v) :
  Base( Rep(v.dimension()) )
/*{\Mcreate introduces the translation by vector $v$.}*/ 
{ int d = v.dimension();
  for (int i = 0; i < d; ++i) {
    ptr()->M_(i,i) = v.homogeneous(d);
    ptr()->M_(i,d) = v.homogeneous(i);
  }
  ptr()->M_(d,d) = v.homogeneous(d);
}

Aff_transformationHd(int d, Scaling, const RT& num, const RT& den) 
  : Base( Rep(d) ) 
/*{\Mcreate returns a scaling by a scale factor |num/den|.}*/
{ Matrix& M = ptr()->M_;
  for (int i = 0; i < d; ++i) M(i,i) = num;
  M(d,d) = den;
}

Aff_transformationHd(int d, Rotation,  
  const RT& sin_num, const RT& cos_num, const RT& den, 
  int e1 = 0, int e2 = 1) : Base( Rep(d) ) 
/*{\Mcreate returns a planar rotation with sine and cosine values
|sin_num/den| and |cos_num/den| in the plane spanned by
the base vectors $b_{e1}$ and $b_{e2}$ in $d$-space. Thus
the default use delivers a planar rotation in the $x$-$y$
plane. \precond $|sin_num|^2 + |cos_num|^2 = |den|^2$
and $0 \leq e_1 < e_2 < d$}*/
{
  CGAL_assertion_msg((sin_num*sin_num + cos_num*cos_num == den*den),
    "planar_rotation: rotation parameters disobey precondition.");
  CGAL_assertion_msg((0<=e1 && e1<=e2 && e2<d),
    "planar_rotation: base vector indices wrong.");
  Matrix& M = ptr()->M_;
  for (int i=0; i<d; i++) M(i,i) = 1;
  M(e1,e1) = cos_num; M(e1,e2) = -sin_num;
  M(e2,e1) = sin_num; M(e2,e2) = cos_num;
  M(d,d) = den;
}

Aff_transformationHd(int d, Rotation, const DirectionHd<RT,LA>& dir, 
  const RT& eps_num, const RT& eps_den, int e1 = 0, int e2 = 1)
/*{\Mcreate returns a planar rotation within the plane spanned by
the base vectors $b_{e1}$ and $b_{e2}$ in $d$-space.  The rotation
parameters are given by the $2$-dimensional direction |dir|, such that
the difference between the sines and cosines of the rotation given by
|dir| and the approximated rotation are at most |num/den| each.\\
\precond |dir.dimension()==2|, |!dir.is_degenerate()| and |num < den|
is positive and $0 \leq e_1 < e_2 < d$ }*/
  : Base( Rep(d) )  
{
  CGAL_assertion(dir.dimension()==2);
  Matrix& M = ptr()->M_;
  for (int i=0; i<d; i++) M(i,i) = RT(1);
  RT sin_num, cos_num, denom;
  rational_rotation_approximation(dir.dx(), dir.dy(),
                                  sin_num, cos_num, denom,
                                  eps_num, eps_den);

  M(e1,e1) = cos_num; M(e1,e2) = -sin_num;
  M(e2,e1) = sin_num; M(e2,e2) = cos_num;
  M(d,d) = denom;
}

/*{\Moperations 5 3}*/

int dimension() const 
{ return ptr()->M_.row_dimension()-1; }
/*{\Mop the dimension of the underlying space }*/

const Matrix& matrix() const { return ptr()->M_; }
/*{\Mop returns the transformation matrix }*/

Vector operator()(const Vector& iv) const
// transforms the ivector by a matrix multiplication
{ return matrix()*iv; }

bool is_odd() const
/*{\Mop returns true iff |\Mvar| is odd.}*/
{ return LA::sign_of_determinant(matrix())<0; }

Aff_transformationHd<RT,LA> inverse() const
/*{\Mop returns the inverse transformation.
\precond |\Mvar.matrix()| is invertible.}*/
{ Aff_transformationHd<RT,LA> Inv; RT D; 
  Vector dummy;
  if ( !LA::inverse(matrix(),Inv.ptr()->M_,D,dummy) ) 
  CGAL_error_msg("Aff_transformationHd::inverse: not invertible.");
  if ( D < 0 ) Inv.ptr()->M_ = -Inv.ptr()->M_;
  return Inv;
}
  
Aff_transformationHd<RT,LA>  
operator*(const Aff_transformationHd<RT,LA>& s) const
/*{\Mbinop composition of transformations. Note that transformations
are not necessarily commutative. |t*s| is the transformation
which transforms first by |t| and then by |s|.}*/
{ CGAL_assertion_msg((dimension()==s.dimension()),
  "Aff_transformationHd::operator*: dimensions disagree.");
  return Aff_transformationHd<RT,LA>(matrix()*s.matrix()); 
}

bool operator==(const Aff_transformationHd<RT,LA>& a1) const
{ if ( this->identical(a1) ) return true;
  return ( matrix() == a1.matrix() );
}
bool operator!=(const Aff_transformationHd<RT,LA>& a1) const
{ return !operator==(a1); }

}; // Aff_transformationHd

template <class RT, class LA>
std::ostream& operator<<(
  std::ostream& os, const Aff_transformationHd<RT,LA>& t) 
{ os << t.matrix(); return os; }

template <class RT, class LA>
std::istream& operator>>(
  std::istream& is, Aff_transformationHd<RT,LA>& t)
{ typename LA::Matrix M(t.dimension());
  is >> M; t = Aff_transformationHd<RT,LA>(M); 
  return is;
}

/*{\Mimplementation 
Affine Transformations are implemented by matrices of integers as an
item type.  All operations like creation, initialization, input and
output on a transformation $t$ take time $O(|t.dimension()|^2)$. |dimension()|
takes constant time.  The operations for inversion and composition
have the cubic costs of the used matrix operations. The space
requirement is $O(|t.dimension()|^2)$. }*/

// ----------------------------- end of file ----------------------------


} //namespace CGAL
#endif // CGAL_AFF_TRANSFORMATIONHD_H
