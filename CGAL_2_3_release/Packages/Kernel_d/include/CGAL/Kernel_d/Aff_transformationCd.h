// ======================================================================
//
// Copyright (c) 2000,2001 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/Kernel_d/Aff_transformationCd.h
// package       : Kernel_d
// maintainer    : Michael Seel <Michael.Seel@mpi-sb.mpg.de>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Michael Seel
// coordinator   : MPI Saarbruecken (Susan.Hert@mpi-sb.mpg.de)
//
// ======================================================================
#ifndef CGAL_AFF_TRANSFORMATIONCD_H
#define CGAL_AFF_TRANSFORMATIONCD_H

#ifndef NOCGALINCL
#include <CGAL/basic.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Handle_for.h>
#include <CGAL/rational_rotation.h>
#endif

CGAL_BEGIN_NAMESPACE

template <class FT, class LA > class Aff_transformationCd;
template <class FT, class LA > class Aff_transformationCd_rep;

template <class FT, class LA>
class Aff_transformationCd_rep : public Ref_counted {
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

public: 
typedef _FT RT;
typedef _FT FT;
typedef _LA LA;
typedef typename _LA::Matrix Matrix;
typedef typename _LA::Vector Vector;

Aff_transformationCd(int d = 0) : Base( Rep(d) ) {}

Aff_transformationCd(int d, Identity_transformation) : Base( Rep(d) )
{ for (int i = 0; i <= d; ++i) ptr->M_(i,i) = FT(1); }

Aff_transformationCd(const Matrix& M) : Base( Rep(M) )
{ CGAL_assertion_msg((M.row_dimension()==M.column_dimension()),
  "Aff_transformationCd:: initialization matrix not quadratic.");
  int d = M.row_dimension(),i;
  for (i=0; i<d-1; ++i) CGAL_assertion(M(d-1,i)==FT(0));
  CGAL_assertion(M(d-1,d-1)==FT(1));
}

#ifndef CGAL_SIMPLE_INTERFACE

template <typename Forward_iterator>
Aff_transformationCd(Scaling, Forward_iterator start, Forward_iterator end) :
  Base( Rep(std::distance(start,end)-1) )
/*{\Mcreate introduces the transformation of $d$-space specified by a
diagonal matrix with entries |set [start,end)| on the diagonal 
(a scaling of the space). \precond |set [start,end)| is a vector of 
dimension $d+1$.}*/
{ int i=0; while (start != end) { ptr->M_(i,i) = *start++;++i; } }

#else
#define FIXATCD(I) \
Aff_transformationCd(Scaling, I start, I end) : \
  Base( Rep(end-start-1) ) \
{ int i=0; while (start != end) { ptr->M_(i,i) = *start++;++i; } }

FIXATCD(int*)
FIXATCD(const int*)
FIXATCD(RT*)
FIXATCD(const RT*)
#undef FIXATCD
#endif

Aff_transformationCd(Translation, const VectorCd<RT,LA>& v) :
  Base( Rep(v.dimension()) )
{ register int d = v.dimension();
  for (int i = 0; i < d; ++i) {
    ptr->M_(i,i) = FT(1);
    ptr->M_(i,d) = v.cartesian(i);
  }
  ptr->M_(d,d) = FT(1);
}

Aff_transformationCd(int d, Scaling, const RT& num, const RT& den) 
  : Base( Rep(d) ) 
{ Matrix& M = ptr->M_;
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
  Matrix& M = ptr->M_;
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
  Matrix& M = ptr->M_;
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
{ return ptr->M_.row_dimension()-1; }

const Matrix& matrix() const { return ptr->M_; }

bool is_odd() const 
{ return LA::sign_of_determinant(matrix())<0; }

Vector operator()(const Vector& v) const
{ CGAL_assertion(matrix().row_dimension()-1==v.dimension());
  Matrix& M = ptr->M_;
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
  Matrix& M = ptr->M_;
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
  if ( !LA::inverse(matrix(),Inv.ptr->M_,D,dummy) ) 
  CGAL_assertion_msg(0,"Aff_transformationCd::inverse: not invertible."); 
  if ( D < FT(0) ) Inv.ptr->M_ = -Inv.ptr->M_;
  return Inv;
}
  
Aff_transformationCd<RT,LA>  
operator*(const Aff_transformationCd<RT,LA>& s) const
{ CGAL_assertion_msg((dimension()==s.dimension()),
  "Aff_transformationCd::operator*: dimensions disagree.");
  return Aff_transformationCd<RT,LA>(matrix()*s.matrix()); 
}

bool operator==(const Aff_transformationCd<RT,LA>& a1) const
{ if ( identical(a1) ) return true;
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


CGAL_END_NAMESPACE
#endif // CGAL_AFF_TRANSFORMATIONCD_H

