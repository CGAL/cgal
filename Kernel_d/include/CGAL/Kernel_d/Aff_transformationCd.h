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
// Author(s)     : Michael Seel

#ifndef CGAL_AFF_TRANSFORMATIONCD_H
#define CGAL_AFF_TRANSFORMATIONCD_H

#include <CGAL/basic.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Handle_for.h>
#include <CGAL/rational_rotation.h>

namespace CGAL {

template <class FT, class LA > class Aff_transformationCd;
template <class FT, class LA > class Aff_transformationCd_rep;

template <class FT, class LA>
class Aff_transformationCd_rep 
{
  friend class Aff_transformationCd<FT,LA>;
  typedef typename LA::Matrix Matrix;
  Matrix M_;
public:
  Aff_transformationCd_rep(int d) : M_(d+1) {}
  Aff_transformationCd_rep(const Matrix& M_init) : M_(M_init) {}
  ~Aff_transformationCd_rep() {}
};

template <class _FT, class _LA>
class Aff_transformationCd : 
  public Handle_for< Aff_transformationCd_rep<_FT,_LA> > { 

  typedef Aff_transformationCd_rep<_FT,_LA> Rep;
  typedef Handle_for<Rep> Base;
  typedef Aff_transformationCd<_FT,_LA> Self;

  using Base::ptr;

public: 
typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename _LA::Matrix Matrix;
typedef typename _LA::Vector Vector;

Aff_transformationCd(int d = 0) : Base( Rep(d) ) {}

Aff_transformationCd(int d, Identity_transformation) : Base( Rep(d) )
{ for (int i = 0; i <= d; ++i) ptr()->M_(i,i) = FT(1); }

Aff_transformationCd(const Matrix& M) : Base( Rep(M) )
{ CGAL_assertion_msg((M.row_dimension()==M.column_dimension()),
  "Aff_transformationCd:: initialization matrix not quadratic.");
  int d = M.row_dimension(),i;
  for (i=0; i<d-1; ++i) CGAL_assertion(M(d-1,i)==FT(0));
  CGAL_assertion(M(d-1,d-1)==FT(1));
}

template <typename Forward_iterator>
Aff_transformationCd(Scaling, Forward_iterator start, Forward_iterator end) :
  Base( Rep(static_cast<int>(std::distance(start,end))-1) )
/*{\Mcreate introduces the transformation of $d$-space specified by a
diagonal matrix with entries |set [start,end)| on the diagonal 
(a scaling of the space). \precond |set [start,end)| is a vector of 
dimension $d+1$.}*/
{ int i=0; while (start != end) { ptr()->M_(i,i) = *start++;++i; } }


Aff_transformationCd(Translation, const VectorCd<RT,LA>& v) :
  Base( Rep(v.dimension()) )
{ int d = v.dimension();
  for (int i = 0; i < d; ++i) {
    ptr()->M_(i,i) = FT(1);
    ptr()->M_(i,d) = v.cartesian(i);
  }
  ptr()->M_(d,d) = FT(1);
}

Aff_transformationCd(int d, Scaling, const RT& num, const RT& den) 
  : Base( Rep(d) ) 
{ Matrix& M = ptr()->M_;
  for (int i = 0; i < d; ++i) M(i,i) = num/den;
  M(d,d) = FT(1);
}

Aff_transformationCd(int d, Rotation,  
  const RT& sin_num, const RT& cos_num, const RT& den, 
  int e1 = 0, int e2 = 1) : Base( Rep(d) )
{
  CGAL_assertion_msg((sin_num*sin_num + cos_num*cos_num == den*den),
    "planar_rotation: rotation parameters disobey precondition.");
  CGAL_assertion_msg((0<=e1 && e1<=e2 && e2<d), 
    "planar_rotation: base vector indices wrong.");
  Matrix& M = ptr()->M_;
  for (int i=0; i<d; i++) M(i,i) = 1;
  M(e1,e1) = cos_num/den; M(e1,e2) = -sin_num/den;
  M(e2,e1) = sin_num/den; M(e2,e2) = cos_num/den;
  M(d,d) = FT(1);
}

Aff_transformationCd(int d, Rotation, const DirectionCd<RT,LA>& dir,
  const RT& eps_num, const RT& eps_den, int e1 = 0, int e2 = 1) 
  : Base( Rep(d) )
{
  CGAL_assertion(dir.dimension()==2);
  Matrix& M = ptr()->M_;
  for (int i=0; i<d; i++) M(i,i) = FT(1);
  RT sin_num, cos_num, denom;
  rational_rotation_approximation(dir.dx(), dir.dy(),
                                  sin_num, cos_num, denom,
                                  eps_num, eps_den);

  M(e1,e1) = cos_num/denom; M(e1,e2) = -sin_num/denom;
  M(e2,e1) = sin_num/denom; M(e2,e2) = cos_num/denom;
  M(d,d) = FT(1);
}

int dimension() const 
{ return ptr()->M_.row_dimension()-1; }

const Matrix& matrix() const { return ptr()->M_; }

bool is_odd() const 
{ return LA::sign_of_determinant(matrix())<0; }

Vector operator()(const Vector& v) const
{ CGAL_assertion(matrix().row_dimension()-1==v.dimension());
  const Matrix& M = ptr()->M_;
  int i,j,d(v.dimension());
  Vector res(d);
  for (i=0; i<d; ++i) { // all rows
    FT cres(0); 
    for (j=0; j<d; ++j) cres+=M(i,j)*v[j]; // per row
    cres += M(i,d);
    res[i]=cres;
  }
  return res; 
}

Vector transform_linearly(const Vector& v) const
{ CGAL_assertion(matrix().row_dimension()-1==v.dimension());
  const Matrix& M = ptr()->M_;
  int i,j,d(v.dimension());
  Vector res(d);
  for (i=0; i<d; ++i) { // all rows
    FT cres(0); 
    for (j=0; j<d; ++j) cres+=M(i,j)*v[j]; // per row
    res[i]=cres;
  }
  return res; 
}


Aff_transformationCd<RT,LA> inverse() const
{ Aff_transformationCd<RT,LA> Inv; RT D; 
  Vector dummy;
  if ( !LA::inverse(matrix(),Inv.ptr()->M_,D,dummy) ) 
  CGAL_error_msg("Aff_transformationCd::inverse: not invertible."); 
  if ( D < FT(0) ) Inv.ptr()->M_ = -Inv.ptr()->M_;
  return Inv;
}
  
Aff_transformationCd<RT,LA>  
operator*(const Aff_transformationCd<RT,LA>& s) const
{ CGAL_assertion_msg((dimension()==s.dimension()),
  "Aff_transformationCd::operator*: dimensions disagree.");
  return Aff_transformationCd<RT,LA>(matrix()*s.matrix()); 
}

bool operator==(const Aff_transformationCd<RT,LA>& a1) const
{ if ( this->identical(a1) ) return true;
  return ( matrix() == a1.matrix() );
}
bool operator!=(const Aff_transformationCd<RT,LA>& a1) const
{ return !operator==(a1); }

}; // Aff_transformationCd

template <class FT, class LA>
std::ostream& operator<<(
  std::ostream& os, const Aff_transformationCd<FT,LA>& t) 
{ os << t.matrix(); return os; }

template <class FT, class LA>
std::istream& operator>>(
  std::istream& is, Aff_transformationCd<FT,LA>& t)
{ typename LA::Matrix M(t.dimension());
  is >> M; t = Aff_transformationCd<FT,LA>(M); 
  return is;
}


} //namespace CGAL
#endif // CGAL_AFF_TRANSFORMATIONCD_H
