#ifndef CGAL_AFF_TRANSFORMATIONCD_H
#define CGAL_AFF_TRANSFORMATIONCD_H

#ifndef NOCGALINCL
#include <CGAL/basic.h>
#include <CGAL/Handle_for.h>
#endif

//-------------------------------------------------------------------------
//
//  NOT FINISHED, JUST CONVERTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
//
//-------------------------------------------------------------------------

CGAL_BEGIN_NAMESPACE

template <class FT, class LA > class Aff_transformationCd;
template <class FT, class LA > class Aff_transformationCd_rep;

template <class FT, class LA>
class Aff_transformationCd_rep : 
  private LA::Matrix, public Ref_counted {
  typedef typename LA::Matrix Matrix;
  friend class Aff_transformationCd<FT,LA>;
public:
  Aff_transformationCd_rep(int d) : Matrix(d+1) {}
  Aff_transformationCd_rep(const Matrix& M_init) : Matrix(M_init) {}
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

Aff_transformationCd(int d = 0, bool identity=false) : Base( Rep(d) )
{ if (identity) for (int i = 0; i <= d; ++i) (*ptr)(i,i) = FT(1); }

Aff_transformationCd(const typename LA::Matrix& M) : Base( Rep(M) )
{ CGAL_assertion_msg((M.row_dimension()==M.column_dimension()),
  "Aff_transformationCd::\
   construction: initialization matrix is not quadratic.");
}

template <typename Forward_iterator>
Aff_transformationCd(Forward_iterator start, Forward_iterator end) :
  Base( Rep(std::distance(start,end)-1) )
/*{\Mcreate introduces the transformation of $d$-space specified by a
diagonal matrix with entries |set [start,end)| on the diagonal 
(a scaling of the space). \precond |set [start,end)| is a vector of 
dimension $d+1$.}*/
{ int i=0; while (start != end) { (*ptr)(i,i) = *start++;++i; } }

Aff_transformationCd(const VectorCd<FT,LA>& v) :
  Base( Rep(v.dimension()) )
/*{\Mcreate introduces the translation by vector $v$.}*/ 
{ register int d = v.dimension();
  for (int i = 0; i < d; ++i) {
    (*ptr)(i,i) = v.homogeneous(d);
    (*ptr)(i,d) = v.homogeneous(i);
  }
  (*ptr)(d,d) = v.homogeneous(d);
}

Aff_transformationCd(int d, const FT& num, const FT& den) : Base( Rep(d) ) 
/*{\Mcreate returns a scaling by a scale factor $\mathit{num}/
\mathit{den}$.}*/
{ typename LA::Matrix& M(*ptr);
  for (int i = 0; i < d; ++i) M(i,i) = num;  
  M(d,d) = den;
}

Aff_transformationCd(int d, 
  const FT& sin_num, const FT& cos_num, const FT& den, 
  int e1 = 0, int e2 = 1); 
/*{\Mcreate returns a planar rotation with sine and cosine values
|sin_num/den| and |cos_num/den| in the plane spanned by
the base vectors $b_{e1}$ and $b_{e2}$ in $d$-space. Thus
the default use delivers a planar rotation in the $x$-$y$
plane. \precond $|sin_num|^2 + |cos_num|^2 = |den|^2$
and $0 \leq e_1 < e_2 < d$}*/

Aff_transformationCd(int d, const DirectionCd<FT,LA>& dir, 
  const FT& num, const FT& den, int e1 = 0, int e2 = 1); 
/*{\Mcreate returns a planar rotation within the plane spanned by
the base vectors $b_{e1}$ and $b_{e2}$ in $d$-space.  The rotation
parameters are given by the $2$-dimensional direction |dir|, such that
the difference between the sines and cosines of the rotation given by
|dir| and the approximated rotation are at most |num/den| each.\\
\precond |dir.dimension()==2|, |!dir.is_degenerate()| and |num < den|
is positive and $0 \leq e_1 < e_2 < d$ }*/

/*{\Moperations 5 3}*/

int dimension() const 
{ return (*ptr).row_dimension()-1; }
/*{\Mop the dimension of the underlying space }*/

const typename LA::Matrix& matrix() const 
{ return static_cast<const typename LA::Matrix&>(*ptr); }
/*{\Mop returns the transformation matrix }*/

typename LA::Vector operator()(const typename LA::Vector& iv) const
// transforms the ivector by a matrix multiplication
{ return matrix()*iv; }


Aff_transformationCd<FT,LA> inverse() const
/*{\Mop returns the inverse transformation.
\precond |\Mvar.matrix()| is invertible.}*/
{ Aff_transformationCd<FT,LA> Inv; FT D; 
  typename LA::Vector dummy;
  if (!LA::inverse((*ptr),(*Inv.ptr),D,dummy)) 
  CGAL_assertion_msg(0,"Aff_transformationCd::inverse: not invertible."); 
  return Inv;
}
  
Aff_transformationCd<FT,LA>  
operator*(const Aff_transformationCd<FT,LA>& s) const
/*{\Mbinop composition of transformations. Note that transformations
are not necessarily commutative. |t*s| is the transformation
which transforms first by |t| and then by |s|.}*/
{ CGAL_assertion_msg((dimension()==s.dimension()),
  "Aff_transformationCd::operator*: dimensions disagree.");
  return Aff_transformationCd<FT,LA>(matrix()*s.matrix()); 
}

bool operator==(const Aff_transformationCd<FT,LA>& a1) const
{ if (identical(a1)) return true;
  return (matrix() == a1.matrix());
}
bool operator!=(const Aff_transformationCd<FT,LA>& a1) const
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

