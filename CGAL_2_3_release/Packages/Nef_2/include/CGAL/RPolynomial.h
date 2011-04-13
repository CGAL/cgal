// ============================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision$
// release_date  : $CGAL_Date$
//
// file          : include/CGAL/RPolynomial.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// source        : nef_2d/RPolynomial.lw
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Polynomials in one variable
// ============================================================================

#ifndef CGAL_RPOLYNOMIAL_H
#define CGAL_RPOLYNOMIAL_H

#include <CGAL/basic.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/IO/io.h>
#undef _DEBUG
#define _DEBUG 3
#include <CGAL/Nef_2/debug.h>

#if defined(_MSC_VER) || defined(__BORLANDC__)
#include <CGAL/Nef_2/vector_MSC.h>
#define CGAL_SIMPLE_NEF_INTERFACE
#define SNIHACK ,char,char
#define SNIINST ,'c','c'
#else
#include <vector>
#define SNIHACK
#define SNIINST
#endif     

class ring_or_field_dont_know {};
class ring_with_gcd {};
class field_with_div {};

template <typename NT>
struct ring_or_field {
  typedef ring_or_field_dont_know kind;
};

template <>
struct ring_or_field<int> {
  typedef ring_with_gcd kind;
  typedef int RT;
  static RT gcd(const RT& a, const RT& b) 
  { if (a == 0)
      if (b == 0)  return 1;
      else         return CGAL_NTS abs(b);
    if (b == 0)    return CGAL_NTS abs(a);
    // here both a and b are non-zero
    int u = CGAL_NTS abs(a);
    int v = CGAL_NTS abs(b);
    if (u < v) v = v%u;
    while (v != 0)
    { int tmp = u % v; 
      u = v;
      v = tmp;
    }
    return u;
  }
};

template <>
struct ring_or_field<long> {
  typedef ring_with_gcd kind;
  typedef long RT;
  static RT gcd(const RT& a, const RT& b) 
  { if (a == 0)
      if (b == 0)  return 1;
      else         return CGAL_NTS abs(b);
    if (b == 0)    return CGAL_NTS abs(a);
    // here both a and b are non-zero
    int u = CGAL_NTS abs(a);
    int v = CGAL_NTS abs(b);
    if (u < v) v = v%u;
    while (v != 0)
    { int tmp = u % v; 
      u = v;
      v = tmp;
    }
    return u;
  }
};


template <>
struct ring_or_field<double> {
  typedef field_with_div kind;
  typedef double RT;
  static RT gcd(const RT&, const RT&) 
  { return 1.0; }
};


CGAL_BEGIN_NAMESPACE

template <class NT> class RPolynomial_rep;

// SPECIALIZE_CLASS(NT,int double) START
// CLASS TEMPLATE NT: 
template <class NT> class RPolynomial;
// CLASS TEMPLATE int: 
CGAL_TEMPLATE_NULL class RPolynomial<int> ;
// CLASS TEMPLATE double: 
CGAL_TEMPLATE_NULL class RPolynomial<double> ;
// SPECIALIZE_CLASS(NT,int double) END

/*{\Mtext \headerline{Range template}}*/
#ifndef CGAL_SIMPLE_NEF_INTERFACE

template <class Forward_iterator>
typename std::iterator_traits<Forward_iterator>::value_type 
gcd_of_range(Forward_iterator its, Forward_iterator ite)
/*{\Mfunc calculates the greates common divisor of the
set of numbers $\{ |*its|, |*++its|, \ldots, |*it| \}$ of type |NT|,
where |++it == ite| and |NT| is the value type of |Forward_iterator|. 
\precond there exists a pairwise gcd operation |NT gcd(NT,NT)| and 
|its!=ite|.}*/
{ CGAL_assertion(its!=ite);
  typedef typename std::iterator_traits<Forward_iterator>::value_type NT;
  NT res = *its++;
  for(; its!=ite; ++its) res = 
    (*its==0 ? res : ring_or_field<NT>::gcd(res, *its));
  if (res==0) res = 1;
  return res;
}

#endif //CGAL_SIMPLE_NEF_INTERFACE


template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
  operator - (const RPolynomial<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  RPolynomial<NT>
  operator + (const RPolynomial<NT>&, const RPolynomial<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  RPolynomial<NT>
  operator - (const RPolynomial<NT>&, const RPolynomial<NT>&);
template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
  operator * (const RPolynomial<NT>&, const RPolynomial<NT>&);
template <class NT> inline RPolynomial<NT>
  operator / (const RPolynomial<NT>&, const RPolynomial<NT>&);

#ifndef _MSC_VER 
template<class NT> /*CGAL_KERNEL_INLINE*/ CGAL::Sign 
  sign(const RPolynomial<NT>& p);
#endif // collides with global CGAL sign

template <class NT> /*CGAL_KERNEL_INLINE*/ double 
  to_double(const RPolynomial<NT>& p) ;
template <class NT>  /*CGAL_KERNEL_INLINE*/ bool 
  is_valid(const RPolynomial<NT>& p) ;
template <class NT> /*CGAL_KERNEL_INLINE*/ bool 
  is_finite(const RPolynomial<NT>& p) ;

template<class NT>  
  std::ostream& operator << (std::ostream& os, const RPolynomial<NT>& p);
template <class NT>  
  std::istream& operator >> (std::istream& is, RPolynomial<NT>& p);

#ifndef _MSC_VER
  // lefthand side
template<class NT> inline RPolynomial<NT> operator + 
  (const NT& num, const RPolynomial<NT>& p2);
template<class NT> inline RPolynomial<NT> operator - 
  (const NT& num, const RPolynomial<NT>& p2);
template<class NT> inline RPolynomial<NT> operator * 
  (const NT& num, const RPolynomial<NT>& p2);
template<class NT> inline RPolynomial<NT> operator / 
  (const NT& num, const RPolynomial<NT>& p2);

  // righthand side
template<class NT> inline RPolynomial<NT> operator + 
  (const RPolynomial<NT>& p1, const NT& num);
template<class NT> inline RPolynomial<NT> operator - 
  (const RPolynomial<NT>& p1, const NT& num);
template<class NT> inline RPolynomial<NT> operator * 
  (const RPolynomial<NT>& p1, const NT& num);
template<class NT> inline RPolynomial<NT> operator / 
  (const RPolynomial<NT>& p1, const NT& num);

  // lefthand side
template<class NT> inline bool operator ==  
  (const NT& num, const RPolynomial<NT>& p);
template<class NT> inline bool operator != 
  (const NT& num, const RPolynomial<NT>& p);
template<class NT> inline bool operator <  
  (const NT& num, const RPolynomial<NT>& p);
template<class NT> inline bool operator <=  
  (const NT& num, const RPolynomial<NT>& p);
template<class NT> inline bool operator >  
  (const NT& num, const RPolynomial<NT>& p);
template<class NT> inline bool operator >=  
  (const NT& num, const RPolynomial<NT>& p);

  // righthand side
template<class NT> inline bool operator ==
  (const RPolynomial<NT>& p, const NT& num);
template<class NT> inline bool operator !=
  (const RPolynomial<NT>& p, const NT& num);
template<class NT> inline bool operator < 
  (const RPolynomial<NT>& p, const NT& num);
template<class NT> inline bool operator <= 
  (const RPolynomial<NT>& p, const NT& num);
template<class NT> inline bool operator > 
  (const RPolynomial<NT>& p, const NT& num);
template<class NT> inline bool operator >=
  (const RPolynomial<NT>& p, const NT& num);

#endif // _MSC_VER

template <class pNT> class RPolynomial_rep : public Ref_counted
{ 
  typedef pNT NT;
  #ifndef CGAL_SIMPLE_NEF_INTERFACE
  typedef std::vector<NT> Vector;
  #else
  typedef CGAL::vector_MSC<NT> Vector;
  #endif
  typedef typename Vector::size_type      size_type;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;
  Vector coeff;

  RPolynomial_rep() : coeff() {}
  RPolynomial_rep(const NT& n) : coeff(1) { coeff[0]=n; }
  RPolynomial_rep(const NT& n, const NT& m) : coeff(2)
    { coeff[0]=n; coeff[1]=m; }
  RPolynomial_rep(const NT& a, const NT& b, const NT& c) : coeff(3)
    { coeff[0]=a; coeff[1]=b; coeff[2]=c; }
  RPolynomial_rep(size_type s) : coeff(s,NT(0)) {}

  #ifndef CGAL_SIMPLE_NEF_INTERFACE

  template <class Forward_iterator>
  RPolynomial_rep(Forward_iterator first, Forward_iterator last SNIHACK) 
    : coeff(first,last) {}

  #else
  template <class Forward_iterator>
  RPolynomial_rep(Forward_iterator first, Forward_iterator last SNIHACK) 
    : coeff() 
  { while (first!=last) coeff.push_back(*first++); }

  #endif

  void reduce() 
  { while ( coeff.size()>1 && coeff.back()==NT(0) ) coeff.pop_back(); }

  friend class RPolynomial<pNT>;
  friend class RPolynomial<int>;
  friend class RPolynomial<double>;
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS  
         (std::istream&, RPolynomial<NT>&);


};

// SPECIALIZE_CLASS(NT,int double) START
// CLASS TEMPLATE NT: 
/*{\Msubst 
typename iterator_traits<Forward_iterator>::value_type#NT
CGAL_NULL_TMPL_ARGS#
}*/

/*{\Manpage{RPolynomial}{NT}{Polynomials in one variable}{p}}*/

template <class pNT> class RPolynomial : 
  public Handle_for< RPolynomial_rep<pNT> >
{
/*{\Mdefinition An instance |\Mvar| of the data type |\Mname| represents
a polynomial $p = a_0 + a_1 x + \ldots a_d x^d $ from the ring |NT[x]|. 
The data type offers standard ring operations and a sign operation which 
determines the sign for the limit process $x \rightarrow \infty$. 

|NT[x]| becomes a unique factorization domain, if the number type
|NT| is either a field type (1) or a unique factorization domain
(2). In both cases there's a polynomial division operation defined.}*/

  /*{\Mtypes 5}*/
  public:
  typedef pNT NT;
  /*{\Mtypemember the component type representing the coefficients.}*/

  typedef Handle_for< RPolynomial_rep<NT> > Base;
  typedef RPolynomial_rep<NT> Rep;
  typedef typename Rep::Vector    Vector;
  typedef typename Rep::size_type size_type;
  typedef typename Rep::iterator  iterator;

  typedef typename Rep::const_iterator const_iterator;
  /*{\Mtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  protected:
  void reduce() { ptr->reduce(); }
  Vector& coeffs() { return ptr->coeff; }
  const Vector& coeffs() const { return ptr->coeff; }
  RPolynomial(size_type s) : Base( RPolynomial_rep<NT>(s) ) {}
  // creates a polynomial of degree s-1

  static NT R_; // for visualization only

  /*{\Mcreation 3}*/
  public:

  RPolynomial()
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| of undefined
  value. }*/
    : Base( RPolynomial_rep<NT>() ) {}
   
  RPolynomial(const NT& a0)
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| representing
  the constant polynomial $a_0$.}*/
    : Base(RPolynomial_rep<NT>(a0)) { reduce(); }

  RPolynomial(const NT& a0, const NT& a1)
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x$. }*/
    : Base(RPolynomial_rep<NT>(a0,a1)) { reduce(); }

  RPolynomial(const NT& a0, const NT& a1,const NT& a2)
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x + a_2 x^2$. }*/
    : Base(RPolynomial_rep<NT>(a0,a1,a2)) { reduce(); }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE
  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last)
  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<NT>(first,last)) { reduce(); }

  #else
  #define RPOL(I)\
  RPolynomial(I first, I last) : \
  Base(RPolynomial_rep<NT>(first,last SNIINST)) { reduce(); }
  RPOL(const NT*)
  // KILL int START
  RPOL(const int*)
  // KILL int END
  // KILL double START
  RPOL(const double*)
  // KILL double END
  #undef RPOL
  #endif // CGAL_SIMPLE_NEF_INTERFACE

  // KILL double START
  RPolynomial(double n) : Base(RPolynomial_rep<NT>(NT(n))) { reduce(); }
  RPolynomial(double n1, double n2) 
    : Base(RPolynomial_rep<NT>(NT(n1),NT(n2))) { reduce(); }
  // KILL double END

  // KILL int START
  RPolynomial(int n) : Base(RPolynomial_rep<NT>(NT(n))) { reduce(); }
  RPolynomial(int n1, int n2)
    : Base(RPolynomial_rep<NT>(NT(n1),NT(n2))) { reduce(); }
  // KILL int END

  RPolynomial(const RPolynomial<NT>& p) : Base(p) {}

  protected: // accessing coefficients internally:
  NT& coeff(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared() && i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  /*{\Moperations 3 3 }*/
  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Mop a random access iterator pointing to $a_0$.}*/
  const_iterator end()   const { return ptr->coeff.end(); }
  /*{\Mop a random access iterator pointing beyond $a_d$.}*/

  int degree() const 
  { return ptr->coeff.size()-1; } 
  /*{\Mop the degree of the polynomial.}*/

  const NT& operator[](unsigned int i) const 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }
  /*{\Marrop the coefficient $a_i$ of the polynomial.}*/

  const NT& operator[](unsigned int i) 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }

  NT eval_at(const NT& r) const
  /*{\Mop evaluates the polynomial at |r|.}*/
  { CGAL_assertion( degree()>=0 );
    NT res = ptr->coeff[0], x = r;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=r; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Mop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { const NT& leading_coeff = ptr->coeff.back();
    if (leading_coeff < NT(0)) return (CGAL::NEGATIVE);
    if (leading_coeff > NT(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Mop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==NT(0); }

  RPolynomial<NT> abs() const
  /*{\Mop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE

  NT content() const
  /*{\Mop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |NT| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }

  #else // CGAL_SIMPLE_NEF_INTERFACE

  NT content() const
  { CGAL_assertion( degree()>=0 );
    iterator its=ptr->coeff.begin(),ite=ptr->coeff.end();
    NT res = *its++;
    for(; its!=ite; ++its) res = 
      (*its==0 ? res : ring_or_field<NT>::gcd(res, *its));
    if (res==0) res = 1;
    return res;
  }

  #endif

  static void set_R(const NT& R) { R_ = R; }

  /*{\Mtext Additionally |\Mname| offers standard arithmetic ring
  opertions like |+,-,*,+=,-=,*=|. By means of the sign operation we can
  also offer comparison predicates as $<,>,\leq,\geq$. Where $p_1 < p_2$
  holds iff $|sign|(p_1 - p_2) < 0$. This data type is fully compliant
  to the requirements of CGAL number types. \setopdims{3cm}{2cm}}*/

  friend  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
    operator - CGAL_NULL_TMPL_ARGS  (const RPolynomial<NT>&);   
                          
  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
    operator + CGAL_NULL_TMPL_ARGS (const RPolynomial<NT>&, 
                                    const RPolynomial<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
    operator - CGAL_NULL_TMPL_ARGS (const RPolynomial<NT>&, 
                                    const RPolynomial<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
    operator * CGAL_NULL_TMPL_ARGS (const RPolynomial<NT>&, 
                                    const RPolynomial<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<NT>  operator / CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2);
  /*{\Mbinopfunc implements polynomial division of |p1| and |p2|. if
  |p1 = p2 * p3| then |p2| is returned. The result is undefined if |p3|
  does not exist in |NT[x]|.  The correct division algorithm is chosen
  according to a traits class |ring_or_field<NT>| provided by the user.
  If |ring_or_field<NT>::kind == ring_with_gcd| then the division is
  done by \emph{pseudo division} based on a |gcd| operation of |NT|.  If
  |ring_or_field<NT>::kind == field_with_div| then the division is done
  by \emph{euclidean division} based on the division operation of the
  field |NT|.

  \textbf{Note} that |NT=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Mtext \headerline{Non member functions}}*/

  static RPolynomial<NT> gcd
    (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2);
  /*{\Mstatic returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |NT=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |NT| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  static void pseudo_div
    (const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
     RPolynomial<NT>& q, RPolynomial<NT>& r, NT& D);
  /*{\Mstatic implements division with remainder on polynomials of 
  the ring |NT[x]|: $D*f = g*q + r$.  \precond |NT| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |NT|.}*/

  static void euclidean_div 
    (const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
     RPolynomial<NT>& q, RPolynomial<NT>& r);
  /*{\Mstatic implements division with remainder on polynomials of 
  the ring |NT[x]|: $f = g*q + r$.  \precond |NT| is a field, i.e.,
  there exists a division operation on |NT|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<NT>& p);


  RPolynomial<NT>& operator += (const RPolynomial<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<NT>& operator -= (const RPolynomial<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) -= p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(-p1[i++]);
    reduce(); return (*this); }

  RPolynomial<NT>& operator *= (const RPolynomial<NT>& p1)
  { (*this)=(*this)*p1; return (*this); }

  RPolynomial<NT>& operator /= (const RPolynomial<NT>& p1)
  { (*this)=(*this)/p1; return (*this); }


  //------------------------------------------------------------------
  // SPECIALIZE_MEMBERS(int double) START
    
  RPolynomial<NT>& operator += (const NT& num)
  { copy_on_write();
    coeff(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const NT& num)
  { copy_on_write();
    coeff(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const NT& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const NT& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (NT)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const int& 
    
  RPolynomial<NT>& operator += (const int& num)
  { copy_on_write();
    coeff(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const int& num)
  { copy_on_write();
    coeff(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (NT)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const double& 
    
  RPolynomial<NT>& operator += (const double& num)
  { copy_on_write();
    coeff(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const double& num)
  { copy_on_write();
    coeff(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (NT)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<NT>& p, const NT& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<NT> s(size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.coeff(i) = b*p[i-k];
    operator-=(s);
  }

};

/*{\Mimplementation This data type is implemented as an item type
via a smart pointer scheme. The coefficients are stored in a vector of
|NT| entries.  The simple arithmetic operations $+,-$ take time
$O(d*T(NT))$, multiplication is quadratic in maximal degree of the
arguments times $T(NT)$, where $T(NT)$ is the time for a corresponding
operation on two instances of the ring type.}*/


// CLASS TEMPLATE int: 
/*{\Xsubst 
 iterator_traits<Forward_iterator>::value_type#int
CGAL_NULL_TMPL_ARGS#
}*/

/*{\Xanpage{RPolynomial}{int}{Polynomials in one variable}{p}}*/

CGAL_TEMPLATE_NULL class RPolynomial<int> : 
  public Handle_for< RPolynomial_rep<int> >
{
/*{\Xdefinition An instance |\Mvar| of the data type |\Mname| represents
a polynomial $p = a_0 + a_1 x + \ldots a_d x^d $ from the ring |int[x]|. 
The data type offers standard ring operations and a sign operation which 
determines the sign for the limit process $x \rightarrow \infty$. 

|int[x]| becomes a unique factorization domain, if the number type
|int| is either a field type (1) or a unique factorization domain
(2). In both cases there's a polynomial division operation defined.}*/

  /*{\Xtypes 5}*/
  public:
  typedef int NT;
  /*{\Xtypemember the component type representing the coefficients.}*/

  typedef Handle_for< RPolynomial_rep<int> > Base;
  typedef RPolynomial_rep<int> Rep;
  typedef  Rep::Vector    Vector;
  typedef  Rep::size_type size_type;
  typedef  Rep::iterator  iterator;

  typedef  Rep::const_iterator const_iterator;
  /*{\Xtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  protected:
  void reduce() { ptr->reduce(); }
  Vector& coeffs() { return ptr->coeff; }
  const Vector& coeffs() const { return ptr->coeff; }
  RPolynomial(size_type s) : Base( RPolynomial_rep<int>(s) ) {}
  // creates a polynomial of degree s-1

  static int R_; // for visualization only

  /*{\Xcreation 3}*/
  public:

  RPolynomial()
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| of undefined
  value. }*/
    : Base( RPolynomial_rep<int>() ) {}
   
  RPolynomial(const int& a0)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the constant polynomial $a_0$.}*/
    : Base(RPolynomial_rep<int>(a0)) { reduce(); }

  RPolynomial(const int& a0, const int& a1)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x$. }*/
    : Base(RPolynomial_rep<int>(a0,a1)) { reduce(); }

  RPolynomial(const int& a0, const int& a1,const int& a2)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x + a_2 x^2$. }*/
    : Base(RPolynomial_rep<int>(a0,a1,a2)) { reduce(); }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE
  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<int>(first,last)) { reduce(); }

  #else
  #define RPOL(I)\
  RPolynomial(I first, I last) : \
  Base(RPolynomial_rep<int>(first,last SNIINST)) { reduce(); }
  RPOL(const int*)
  // KILL double START
  RPOL(const double*)
  // KILL double END
  #undef RPOL
  #endif // CGAL_SIMPLE_NEF_INTERFACE

  // KILL double START
  RPolynomial(double n) : Base(RPolynomial_rep<int>(int(n))) { reduce(); }
  RPolynomial(double n1, double n2) 
    : Base(RPolynomial_rep<int>(int(n1),int(n2))) { reduce(); }
  // KILL double END

  RPolynomial(const RPolynomial<int>& p) : Base(p) {}

  protected: // accessing coefficients internally:
  int& coeff(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared() && i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  /*{\Xoperations 3 3 }*/
  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Xop a random access iterator pointing to $a_0$.}*/
  const_iterator end()   const { return ptr->coeff.end(); }
  /*{\Xop a random access iterator pointing beyond $a_d$.}*/

  int degree() const 
  { return ptr->coeff.size()-1; } 
  /*{\Xop the degree of the polynomial.}*/

  const int& operator[](unsigned int i) const 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }
  /*{\Xarrop the coefficient $a_i$ of the polynomial.}*/

  const int& operator[](unsigned int i) 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }

  int eval_at(const int& r) const
  /*{\Xop evaluates the polynomial at |r|.}*/
  { CGAL_assertion( degree()>=0 );
    int res = ptr->coeff[0], x = r;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=r; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Xop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { const int& leading_coeff = ptr->coeff.back();
    if (leading_coeff < int(0)) return (CGAL::NEGATIVE);
    if (leading_coeff > int(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Xop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==int(0); }

  RPolynomial<int> abs() const
  /*{\Xop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE

  int content() const
  /*{\Xop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |int| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }

  #else // CGAL_SIMPLE_NEF_INTERFACE

  int content() const
  { CGAL_assertion( degree()>=0 );
    iterator its=ptr->coeff.begin(),ite=ptr->coeff.end();
    int res = *its++;
    for(; its!=ite; ++its) res = 
      (*its==0 ? res : ring_or_field<int>::gcd(res, *its));
    if (res==0) res = 1;
    return res;
  }

  #endif

  static void set_R(const int& R) { R_ = R; }

  /*{\Xtext Additionally |\Mname| offers standard arithmetic ring
  opertions like |+,-,*,+=,-=,*=|. By means of the sign operation we can
  also offer comparison predicates as $<,>,\leq,\geq$. Where $p_1 < p_2$
  holds iff $|sign|(p_1 - p_2) < 0$. This data type is fully compliant
  to the requirements of CGAL number types. \setopdims{3cm}{2cm}}*/

  friend  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<int>
    operator - CGAL_NULL_TMPL_ARGS  (const RPolynomial<int>&);   
                          
  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<int>
    operator + CGAL_NULL_TMPL_ARGS (const RPolynomial<int>&, 
                                    const RPolynomial<int>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<int>
    operator - CGAL_NULL_TMPL_ARGS (const RPolynomial<int>&, 
                                    const RPolynomial<int>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<int>
    operator * CGAL_NULL_TMPL_ARGS (const RPolynomial<int>&, 
                                    const RPolynomial<int>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<int>  operator / CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<int>& p1, const RPolynomial<int>& p2);
  /*{\Xbinopfunc implements polynomial division of |p1| and |p2|. if
  |p1 = p2 * p3| then |p2| is returned. The result is undefined if |p3|
  does not exist in |int[x]|.  The correct division algorithm is chosen
  according to a traits class |ring_or_field<int>| provided by the user.
  If |ring_or_field<int>::kind == ring_with_gcd| then the division is
  done by \emph{pseudo division} based on a |gcd| operation of |int|.  If
  |ring_or_field<int>::kind == field_with_div| then the division is done
  by \emph{euclidean division} based on the division operation of the
  field |int|.

  \textbf{Note} that |int=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Xtext \headerline{Non member functions}}*/

  static RPolynomial<int> gcd
    (const RPolynomial<int>& p1, const RPolynomial<int>& p2);
  /*{\Xstatic returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |int=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |int| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  static void pseudo_div
    (const RPolynomial<int>& f, const RPolynomial<int>& g, 
     RPolynomial<int>& q, RPolynomial<int>& r, int& D);
  /*{\Xstatic implements division with remainder on polynomials of 
  the ring |int[x]|: $D*f = g*q + r$.  \precond |int| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |int|.}*/

  static void euclidean_div 
    (const RPolynomial<int>& f, const RPolynomial<int>& g, 
     RPolynomial<int>& q, RPolynomial<int>& r);
  /*{\Xstatic implements division with remainder on polynomials of 
  the ring |int[x]|: $f = g*q + r$.  \precond |int| is a field, i.e.,
  there exists a division operation on |int|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<int>& p);


  RPolynomial<int>& operator += (const RPolynomial<int>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<int>& operator -= (const RPolynomial<int>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) -= p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(-p1[i++]);
    reduce(); return (*this); }

  RPolynomial<int>& operator *= (const RPolynomial<int>& p1)
  { (*this)=(*this)*p1; return (*this); }

  RPolynomial<int>& operator /= (const RPolynomial<int>& p1)
  { (*this)=(*this)/p1; return (*this); }


  //------------------------------------------------------------------
  // SPECIALIZE_MEMBERS(int double) START
    
  RPolynomial<int>& operator += (const int& num)
  { copy_on_write();
    coeff(0) += (int)num; return *this; }

  RPolynomial<int>& operator -= (const int& num)
  { copy_on_write();
    coeff(0) -= (int)num; return *this; }

  RPolynomial<int>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (int)num; 
    reduce(); return *this; }

  RPolynomial<int>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (int)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const double& 
    
  RPolynomial<int>& operator += (const double& num)
  { copy_on_write();
    coeff(0) += (int)num; return *this; }

  RPolynomial<int>& operator -= (const double& num)
  { copy_on_write();
    coeff(0) -= (int)num; return *this; }

  RPolynomial<int>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (int)num; 
    reduce(); return *this; }

  RPolynomial<int>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (int)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<int>& p, const int& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<int> s(size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.coeff(i) = b*p[i-k];
    operator-=(s);
  }

};

/*{\Ximplementation This data type is implemented as an item type
via a smart pointer scheme. The coefficients are stored in a vector of
|int| entries.  The simple arithmetic operations $+,-$ take time
$O(d*T(int))$, multiplication is quadratic in maximal degree of the
arguments times $T(int)$, where $T(int)$ is the time for a corresponding
operation on two instances of the ring type.}*/


// CLASS TEMPLATE double: 
/*{\Xsubst 
 iterator_traits<Forward_iterator>::value_type#double
CGAL_NULL_TMPL_ARGS#
}*/

/*{\Xanpage{RPolynomial}{double}{Polynomials in one variable}{p}}*/

CGAL_TEMPLATE_NULL class RPolynomial<double> : 
  public Handle_for< RPolynomial_rep<double> >
{
/*{\Xdefinition An instance |\Mvar| of the data type |\Mname| represents
a polynomial $p = a_0 + a_1 x + \ldots a_d x^d $ from the ring |double[x]|. 
The data type offers standard ring operations and a sign operation which 
determines the sign for the limit process $x \rightarrow \infty$. 

|double[x]| becomes a unique factorization domain, if the number type
|double| is either a field type (1) or a unique factorization domain
(2). In both cases there's a polynomial division operation defined.}*/

  /*{\Xtypes 5}*/
  public:
  typedef double NT;
  /*{\Xtypemember the component type representing the coefficients.}*/

  typedef Handle_for< RPolynomial_rep<double> > Base;
  typedef RPolynomial_rep<double> Rep;
  typedef  Rep::Vector    Vector;
  typedef  Rep::size_type size_type;
  typedef  Rep::iterator  iterator;

  typedef  Rep::const_iterator const_iterator;
  /*{\Xtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  protected:
  void reduce() { ptr->reduce(); }
  Vector& coeffs() { return ptr->coeff; }
  const Vector& coeffs() const { return ptr->coeff; }
  RPolynomial(size_type s) : Base( RPolynomial_rep<double>(s) ) {}
  // creates a polynomial of degree s-1

  static double R_; // for visualization only

  /*{\Xcreation 3}*/
  public:

  RPolynomial()
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| of undefined
  value. }*/
    : Base( RPolynomial_rep<double>() ) {}
   
  RPolynomial(const double& a0)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the constant polynomial $a_0$.}*/
    : Base(RPolynomial_rep<double>(a0)) { reduce(); }

  RPolynomial(const double& a0, const double& a1)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x$. }*/
    : Base(RPolynomial_rep<double>(a0,a1)) { reduce(); }

  RPolynomial(const double& a0, const double& a1,const double& a2)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial $a_0 + a_1 x + a_2 x^2$. }*/
    : Base(RPolynomial_rep<double>(a0,a1,a2)) { reduce(); }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE
  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last)
  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<double>(first,last)) { reduce(); }

  #else
  #define RPOL(I)\
  RPolynomial(I first, I last) : \
  Base(RPolynomial_rep<double>(first,last SNIINST)) { reduce(); }
  RPOL(const double*)
  // KILL int START
  RPOL(const int*)
  // KILL int END
  #undef RPOL
  #endif // CGAL_SIMPLE_NEF_INTERFACE

  // KILL int START
  RPolynomial(int n) : Base(RPolynomial_rep<double>(double(n))) { reduce(); }
  RPolynomial(int n1, int n2)
    : Base(RPolynomial_rep<double>(double(n1),double(n2))) { reduce(); }
  // KILL int END

  RPolynomial(const RPolynomial<double>& p) : Base(p) {}

  protected: // accessing coefficients internally:
  double& coeff(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared() && i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  /*{\Xoperations 3 3 }*/
  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Xop a random access iterator pointing to $a_0$.}*/
  const_iterator end()   const { return ptr->coeff.end(); }
  /*{\Xop a random access iterator pointing beyond $a_d$.}*/

  int degree() const 
  { return ptr->coeff.size()-1; } 
  /*{\Xop the degree of the polynomial.}*/

  const double& operator[](unsigned int i) const 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }
  /*{\Xarrop the coefficient $a_i$ of the polynomial.}*/

  const double& operator[](unsigned int i) 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }

  double eval_at(const double& r) const
  /*{\Xop evaluates the polynomial at |r|.}*/
  { CGAL_assertion( degree()>=0 );
    double res = ptr->coeff[0], x = r;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=r; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Xop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { const double& leading_coeff = ptr->coeff.back();
    if (leading_coeff < double(0)) return (CGAL::NEGATIVE);
    if (leading_coeff > double(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Xop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==double(0); }

  RPolynomial<double> abs() const
  /*{\Xop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  #ifndef CGAL_SIMPLE_NEF_INTERFACE

  double content() const
  /*{\Xop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |double| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }

  #else // CGAL_SIMPLE_NEF_INTERFACE

  double content() const
  { CGAL_assertion( degree()>=0 );
    iterator its=ptr->coeff.begin(),ite=ptr->coeff.end();
    double res = *its++;
    for(; its!=ite; ++its) res = 
      (*its==0 ? res : ring_or_field<double>::gcd(res, *its));
    if (res==0) res = 1;
    return res;
  }

  #endif

  static void set_R(const double& R) { R_ = R; }

  /*{\Xtext Additionally |\Mname| offers standard arithmetic ring
  opertions like |+,-,*,+=,-=,*=|. By means of the sign operation we can
  also offer comparison predicates as $<,>,\leq,\geq$. Where $p_1 < p_2$
  holds iff $|sign|(p_1 - p_2) < 0$. This data type is fully compliant
  to the requirements of CGAL number types. \setopdims{3cm}{2cm}}*/

  friend  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<double>
    operator - CGAL_NULL_TMPL_ARGS  (const RPolynomial<double>&);   
                          
  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<double>
    operator + CGAL_NULL_TMPL_ARGS (const RPolynomial<double>&, 
                                    const RPolynomial<double>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<double>
    operator - CGAL_NULL_TMPL_ARGS (const RPolynomial<double>&, 
                                    const RPolynomial<double>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<double>
    operator * CGAL_NULL_TMPL_ARGS (const RPolynomial<double>&, 
                                    const RPolynomial<double>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<double>  operator / CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<double>& p1, const RPolynomial<double>& p2);
  /*{\Xbinopfunc implements polynomial division of |p1| and |p2|. if
  |p1 = p2 * p3| then |p2| is returned. The result is undefined if |p3|
  does not exist in |double[x]|.  The correct division algorithm is chosen
  according to a traits class |ring_or_field<double>| provided by the user.
  If |ring_or_field<double>::kind == ring_with_gcd| then the division is
  done by \emph{pseudo division} based on a |gcd| operation of |double|.  If
  |ring_or_field<double>::kind == field_with_div| then the division is done
  by \emph{euclidean division} based on the division operation of the
  field |double|.

  \textbf{Note} that |double=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Xtext \headerline{Non member functions}}*/

  static RPolynomial<double> gcd
    (const RPolynomial<double>& p1, const RPolynomial<double>& p2);
  /*{\Xstatic returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |double=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |double| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  static void pseudo_div
    (const RPolynomial<double>& f, const RPolynomial<double>& g, 
     RPolynomial<double>& q, RPolynomial<double>& r, double& D);
  /*{\Xstatic implements division with remainder on polynomials of 
  the ring |double[x]|: $D*f = g*q + r$.  \precond |double| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |double|.}*/

  static void euclidean_div 
    (const RPolynomial<double>& f, const RPolynomial<double>& g, 
     RPolynomial<double>& q, RPolynomial<double>& r);
  /*{\Xstatic implements division with remainder on polynomials of 
  the ring |double[x]|: $f = g*q + r$.  \precond |double| is a field, i.e.,
  there exists a division operation on |double|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<double>& p);


  RPolynomial<double>& operator += (const RPolynomial<double>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<double>& operator -= (const RPolynomial<double>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) -= p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(-p1[i++]);
    reduce(); return (*this); }

  RPolynomial<double>& operator *= (const RPolynomial<double>& p1)
  { (*this)=(*this)*p1; return (*this); }

  RPolynomial<double>& operator /= (const RPolynomial<double>& p1)
  { (*this)=(*this)/p1; return (*this); }


  //------------------------------------------------------------------
  // SPECIALIZE_MEMBERS(int double) START
    
  RPolynomial<double>& operator += (const double& num)
  { copy_on_write();
    coeff(0) += (double)num; return *this; }

  RPolynomial<double>& operator -= (const double& num)
  { copy_on_write();
    coeff(0) -= (double)num; return *this; }

  RPolynomial<double>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (double)num; 
    reduce(); return *this; }

  RPolynomial<double>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (double)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const int& 
    
  RPolynomial<double>& operator += (const int& num)
  { copy_on_write();
    coeff(0) += (double)num; return *this; }

  RPolynomial<double>& operator -= (const int& num)
  { copy_on_write();
    coeff(0) -= (double)num; return *this; }

  RPolynomial<double>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= (double)num; 
    reduce(); return *this; }

  RPolynomial<double>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= (double)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<double>& p, const double& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<double> s(size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.coeff(i) = b*p[i-k];
    operator-=(s);
  }

};

/*{\Ximplementation This data type is implemented as an item type
via a smart pointer scheme. The coefficients are stored in a vector of
|double| entries.  The simple arithmetic operations $+,-$ take time
$O(d*T(double))$, multiplication is quadratic in maximal degree of the
arguments times $T(double)$, where $T(double)$ is the time for a corresponding
operation on two instances of the ring type.}*/


// SPECIALIZE_CLASS(NT,int double) END

template <class NT> NT RPolynomial<NT>::R_;
int    RPolynomial<int>::R_;
double RPolynomial<double>::R_;

template <class NT> /*CGAL_KERNEL_INLINE*/ double to_double 
  (const RPolynomial<NT>& p) 
  { return (CGAL::to_double(p.eval_at(RPolynomial<NT>::R_))); }

template <class NT>  /*CGAL_KERNEL_INLINE*/ bool is_valid 
  (const RPolynomial<NT>& p) 
  { return (CGAL::is_valid(p[0])); }

template <class NT> /*CGAL_KERNEL_INLINE*/ bool is_finite 
  (const RPolynomial<NT>& p) 
  { return CGAL::is_finite(p[0]); }

template <class NT> /*CGAL_KERNEL_INLINE*/ CGAL::io_Operator 
  io_tag(const RPolynomial<NT>&) 
  { return CGAL::io_Operator(); }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator - (const RPolynomial<NT>& p)
{
  CGAL_assertion(p.degree()>=0);
  RPolynomial<NT> res(p.coeffs().begin(),p.coeffs().end());
  typename RPolynomial<NT>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator + (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{ 
  typedef typename RPolynomial<NT>::size_type size_type;
  CGAL_assertion(p1.degree()>=0 && p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial<NT>  p( size_type(max + 1));
  for (i = 0; i <= min; ++i ) p.coeff(i) = p1[i]+p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)=p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator - (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{ 
  typedef typename RPolynomial<NT>::size_type size_type;
  CGAL_assertion(p1.degree()>=0 && p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial<NT>  p( size_type(max+1) );
  for (i = 0; i <= min; ++i ) p.coeff(i)=p1[i]-p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)= -p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=  p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator * (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{
  typedef typename RPolynomial<NT>::size_type size_type;
  CGAL_assertion(p1.degree()>=0 && p2.degree()>=0);
  RPolynomial<NT>  p( size_type(p1.degree()+p2.degree()+1) ); 
  // initialized with zeros
  for (int i=0; i <= p1.degree(); ++i)
    for (int j=0; j <= p2.degree(); ++j)
      p.coeff(i+j) += (p1[i]*p2[j]); 
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> divop (const RPolynomial<NT>& p1, 
                       const RPolynomial<NT>& p2,
                       ring_or_field_dont_know)
{
  CGAL_assertion_msg(0,"\n\
  The division operation on polynomials requires that you\n\
  specify if your number type provides a binary gcd() operation\n\
  or is a field type including an operator/().\n\
  You do this by creating a specialized class:\n\
  template <> class ring_or_field<yourNT> with a member type:\n\
  typedef ring_with_gcd kind; OR\n\
  typedef field_with_div kind;\n");
  return RPolynomial<NT>(); // never reached
}


template <class NT> inline
RPolynomial<NT> operator / (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{ return divop(p1,p2,ring_or_field<NT>::kind()); }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> divop (const RPolynomial<NT>& p1, 
                       const RPolynomial<NT>& p2,
                       field_with_div)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  RPolynomial<NT> q,r;
  RPolynomial<NT>::euclidean_div(p1,p2,q,r);
  CGAL_postcondition( (p2*q+r==p1) );
  return q;
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> divop (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2,
                       ring_with_gcd)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  RPolynomial<NT> q,r; NT D; 
  RPolynomial<NT>::pseudo_div(p1,p2,q,r,D); 
  CGAL_postcondition( (p2*q+r==p1*RPolynomial<NT>(D)) );
  return q/=D;
}


template <class NT> 
inline RPolynomial<NT> 
gcd(const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
{ return RPolynomial<NT>::gcd(p1,p2); }

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator == 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator != 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <  
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::NEGATIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <= 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >  
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >= 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::NEGATIVE ); }    

#ifndef _MSC_VER 
template <class NT> /*CGAL_KERNEL_INLINE*/ CGAL::Sign 
  sign(const RPolynomial<NT>& p)
  { return p.sign(); }
#endif // collides with global CGAL sign

#ifndef _MSC_VER
//------------------------------------------------------------------
// SPECIALIZE_FUNCTION(NT,int double) START
// SPECIALIZING inline to :

  // lefthand side
  inline    RPolynomial<int> operator + 
  (const int& num, const RPolynomial<int>& p2)
  { return (RPolynomial<int>(num) + p2); }
  inline    RPolynomial<int> operator - 
  (const int& num, const RPolynomial<int>& p2)
  { return (RPolynomial<int>(num) - p2); }
  inline    RPolynomial<int> operator * 
  (const int& num, const RPolynomial<int>& p2)
  { return (RPolynomial<int>(num) * p2); }
  inline    RPolynomial<int> operator / 
  (const int& num, const RPolynomial<int>& p2)
  { return (RPolynomial<int>(num)/p2); }

  // righthand side
  inline    RPolynomial<int> operator + 
  (const RPolynomial<int>& p1, const int& num)
  { return (p1 + RPolynomial<int>(num)); }
  inline    RPolynomial<int> operator - 
  (const RPolynomial<int>& p1, const int& num)
  { return (p1 - RPolynomial<int>(num)); }
  inline    RPolynomial<int> operator * 
  (const RPolynomial<int>& p1, const int& num)
  { return (p1 * RPolynomial<int>(num)); }
  inline    RPolynomial<int> operator / 
  (const RPolynomial<int>& p1, const int& num)
  { return (p1 / RPolynomial<int>(num)); }


  // lefthand side
  inline    bool operator ==  
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() == CGAL::ZERO );}
  inline    bool operator != 
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() != CGAL::ZERO );}
  inline    bool operator <  
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() == CGAL::NEGATIVE );}
  inline    bool operator <=  
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() != CGAL::POSITIVE );}
  inline    bool operator >  
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() == CGAL::POSITIVE );}
  inline    bool operator >=  
  (const int& num, const RPolynomial<int>& p) 
  { return ( (RPolynomial<int>(num)-p).sign() != CGAL::NEGATIVE );}

  // righthand side
  inline    bool operator ==
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() == CGAL::ZERO );}
  inline    bool operator !=
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() != CGAL::ZERO );}
  inline    bool operator < 
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() == CGAL::NEGATIVE );}
  inline    bool operator <= 
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() != CGAL::POSITIVE );}
  inline    bool operator > 
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() == CGAL::POSITIVE );}
  inline    bool operator >=
  (const RPolynomial<int>& p, const int& num) 
  { return ( (p-RPolynomial<int>(num)).sign() != CGAL::NEGATIVE );}

// SPECIALIZING pure int params:

  // lefthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const int& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) + p2); }
  template <class NT>    RPolynomial<NT> operator - 
  (const int& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) - p2); }
  template <class NT>    RPolynomial<NT> operator * 
  (const int& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) * p2); }
  template <class NT>    RPolynomial<NT> operator / 
  (const int& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num)/p2); }

  // righthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const RPolynomial<NT>& p1, const int& num)
  { return (p1 + RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator - 
  (const RPolynomial<NT>& p1, const int& num)
  { return (p1 - RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator * 
  (const RPolynomial<NT>& p1, const int& num)
  { return (p1 * RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator / 
  (const RPolynomial<NT>& p1, const int& num)
  { return (p1 / RPolynomial<NT>(num)); }


  // lefthand side
  template <class NT>    bool operator ==  
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::ZERO );}
  template <class NT>    bool operator != 
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::ZERO );}
  template <class NT>    bool operator <  
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <=  
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator >  
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=  
  (const int& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::NEGATIVE );}

  // righthand side
  template <class NT>    bool operator ==
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::ZERO );}
  template <class NT>    bool operator !=
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::ZERO );}
  template <class NT>    bool operator < 
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <= 
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator > 
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=
  (const RPolynomial<NT>& p, const int& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::NEGATIVE );}

// SPECIALIZING inline to :

  // lefthand side
  inline    RPolynomial<double> operator + 
  (const double& num, const RPolynomial<double>& p2)
  { return (RPolynomial<double>(num) + p2); }
  inline    RPolynomial<double> operator - 
  (const double& num, const RPolynomial<double>& p2)
  { return (RPolynomial<double>(num) - p2); }
  inline    RPolynomial<double> operator * 
  (const double& num, const RPolynomial<double>& p2)
  { return (RPolynomial<double>(num) * p2); }
  inline    RPolynomial<double> operator / 
  (const double& num, const RPolynomial<double>& p2)
  { return (RPolynomial<double>(num)/p2); }

  // righthand side
  inline    RPolynomial<double> operator + 
  (const RPolynomial<double>& p1, const double& num)
  { return (p1 + RPolynomial<double>(num)); }
  inline    RPolynomial<double> operator - 
  (const RPolynomial<double>& p1, const double& num)
  { return (p1 - RPolynomial<double>(num)); }
  inline    RPolynomial<double> operator * 
  (const RPolynomial<double>& p1, const double& num)
  { return (p1 * RPolynomial<double>(num)); }
  inline    RPolynomial<double> operator / 
  (const RPolynomial<double>& p1, const double& num)
  { return (p1 / RPolynomial<double>(num)); }


  // lefthand side
  inline    bool operator ==  
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() == CGAL::ZERO );}
  inline    bool operator != 
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() != CGAL::ZERO );}
  inline    bool operator <  
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() == CGAL::NEGATIVE );}
  inline    bool operator <=  
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() != CGAL::POSITIVE );}
  inline    bool operator >  
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() == CGAL::POSITIVE );}
  inline    bool operator >=  
  (const double& num, const RPolynomial<double>& p) 
  { return ( (RPolynomial<double>(num)-p).sign() != CGAL::NEGATIVE );}

  // righthand side
  inline    bool operator ==
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() == CGAL::ZERO );}
  inline    bool operator !=
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() != CGAL::ZERO );}
  inline    bool operator < 
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() == CGAL::NEGATIVE );}
  inline    bool operator <= 
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() != CGAL::POSITIVE );}
  inline    bool operator > 
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() == CGAL::POSITIVE );}
  inline    bool operator >=
  (const RPolynomial<double>& p, const double& num) 
  { return ( (p-RPolynomial<double>(num)).sign() != CGAL::NEGATIVE );}

// SPECIALIZING pure double params:

  // lefthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const double& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) + p2); }
  template <class NT>    RPolynomial<NT> operator - 
  (const double& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) - p2); }
  template <class NT>    RPolynomial<NT> operator * 
  (const double& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) * p2); }
  template <class NT>    RPolynomial<NT> operator / 
  (const double& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num)/p2); }

  // righthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const RPolynomial<NT>& p1, const double& num)
  { return (p1 + RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator - 
  (const RPolynomial<NT>& p1, const double& num)
  { return (p1 - RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator * 
  (const RPolynomial<NT>& p1, const double& num)
  { return (p1 * RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator / 
  (const RPolynomial<NT>& p1, const double& num)
  { return (p1 / RPolynomial<NT>(num)); }


  // lefthand side
  template <class NT>    bool operator ==  
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::ZERO );}
  template <class NT>    bool operator != 
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::ZERO );}
  template <class NT>    bool operator <  
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <=  
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator >  
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=  
  (const double& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::NEGATIVE );}

  // righthand side
  template <class NT>    bool operator ==
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::ZERO );}
  template <class NT>    bool operator !=
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::ZERO );}
  template <class NT>    bool operator < 
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <= 
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator > 
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=
  (const RPolynomial<NT>& p, const double& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::NEGATIVE );}

// SPECIALIZE_FUNCTION ORIGINAL

  // lefthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const NT& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) + p2); }
  template <class NT>    RPolynomial<NT> operator - 
  (const NT& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) - p2); }
  template <class NT>    RPolynomial<NT> operator * 
  (const NT& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num) * p2); }
  template <class NT>    RPolynomial<NT> operator / 
  (const NT& num, const RPolynomial<NT>& p2)
  { return (RPolynomial<NT>(num)/p2); }

  // righthand side
  template <class NT>    RPolynomial<NT> operator + 
  (const RPolynomial<NT>& p1, const NT& num)
  { return (p1 + RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator - 
  (const RPolynomial<NT>& p1, const NT& num)
  { return (p1 - RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator * 
  (const RPolynomial<NT>& p1, const NT& num)
  { return (p1 * RPolynomial<NT>(num)); }
  template <class NT>    RPolynomial<NT> operator / 
  (const RPolynomial<NT>& p1, const NT& num)
  { return (p1 / RPolynomial<NT>(num)); }


  // lefthand side
  template <class NT>    bool operator ==  
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::ZERO );}
  template <class NT>    bool operator != 
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::ZERO );}
  template <class NT>    bool operator <  
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <=  
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator >  
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=  
  (const NT& num, const RPolynomial<NT>& p) 
  { return ( (RPolynomial<NT>(num)-p).sign() != CGAL::NEGATIVE );}

  // righthand side
  template <class NT>    bool operator ==
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::ZERO );}
  template <class NT>    bool operator !=
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::ZERO );}
  template <class NT>    bool operator < 
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::NEGATIVE );}
  template <class NT>    bool operator <= 
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::POSITIVE );}
  template <class NT>    bool operator > 
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() == CGAL::POSITIVE );}
  template <class NT>    bool operator >=
  (const RPolynomial<NT>& p, const NT& num) 
  { return ( (p-RPolynomial<NT>(num)).sign() != CGAL::NEGATIVE );}

// SPECIALIZE_FUNCTION(NT,int double) END
//------------------------------------------------------------------
#endif // _MSC_VER CGAL_CFG_MATCHING_BUG_2


template <class NT> 
void print_monomial(std::ostream& os, const NT& n, int i)
{
  if (i==0) os << n;
  if (i==1) os << n << "R";
  if (i>1)  os << n << "R^" << i;
}

#define RPOLYNOMIAL_EXPLICIT_OUTPUT

// I/O 
template <class NT>
std::ostream& operator << (std::ostream& os, const RPolynomial<NT>& p)
{
  int i;
  switch( os.iword(CGAL::IO::mode) )
  {
    case CGAL::IO::ASCII :
      os << p.degree() << ' ';
      for(i=0; i<=p.degree(); ++i) 
        os << p[i] << ' ';
      return os;
    case CGAL::IO::BINARY :
      CGAL::write(os, p.degree());
      for(i=0; i<=p.degree(); ++i) 
        CGAL::write(os, p[i]);
      return os;
    default: 
#ifndef RPOLYNOMIAL_EXPLICIT_OUTPUT
      os << "RPolynomial(" << p.degree() << ", ";
      for(i=0; i<=p.degree(); ++i) {
        os << p[i];
        if (i < p.degree()) os << ", ";
      }
      os << ")";
#else
      print_monomial(os,p[p.degree()],p.degree());
      for(i=p.degree()-1; i>=0; --i) {
        if (p[i]!=NT(0)) { os << " + "; print_monomial(os,p[i],i); }
      }    
#endif
      return os;
  }
}

template <class NT>
std::istream& operator >> (std::istream& is, RPolynomial<NT>& p)
{ 
  int i,d;
  NT c;
  switch( is.iword(CGAL::IO::mode) )
  { 
    case CGAL::IO::ASCII : 
      is >> d;
      if (d < 0) p = RPolynomial<NT>();
      else {
        typename RPolynomial<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) is >> coeffs[i];
        p = RPolynomial<NT>(coeffs.begin(),coeffs.end());
      }
      break;
    case CGAL::IO::BINARY :
      CGAL::read(is, d);
      if (d < 0) p = RPolynomial<NT>();
      else {
        typename RPolynomial<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) 
        { CGAL::read(is,c); coeffs[i]=c; }
        p = RPolynomial<NT>(coeffs.begin(),coeffs.end());
      }
      break;
    default:
      CGAL_assertion_msg(0,"\nStream must be in ascii or binary mode\n");
      break;
  }
  return is;
}



// SPECIALIZE_IMPLEMENTATION(NT,int double) START
// SPECIALIZING to :
 /*CGAL_KERNEL_MEDIUM_INLINE*/ 
void RPolynomial<int>::euclidean_div(
  const RPolynomial<int>& f, const RPolynomial<int>& g,
  RPolynomial<int>& q, RPolynomial<int>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = RPolynomial<int>(int(0)); }
  else { qd = rd-gd+1; q = RPolynomial<int>(size_t(qd)); }
  while ( rd >= gd ) {
    int S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}


 /*CGAL_KERNEL_MEDIUM_INLINE*/   
void RPolynomial<int>::pseudo_div(
  const RPolynomial<int>& f, const RPolynomial<int>& g, 
  RPolynomial<int>& q, RPolynomial<int>& r, int& D)
{
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = RPolynomial<int>(0); r = f; D = 1; 
    CGAL_postcondition(RPolynomial<int>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = RPolynomial<int>( size_t(delta) );
  int G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  RPolynomial<int> res = RPolynomial<int>(D)*f;
  TRACEN("  pseudo_div start "<<res<<" "<<qd<<" "<<q.degree());
  while (qd >= 0) {
    int F = res[rd]; // highest order coeff of res
    int t = F/G;     // ensured to be integer by multiplication of D
    q.coeff(qd) = t;    // store q coeff
    res.minus_offsetmult(g,t,qd); 
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res;
  CGAL_postcondition(RPolynomial<int>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}


 /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<int> RPolynomial<int>::gcd(
  const RPolynomial<int>& p1, const RPolynomial<int>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return RPolynomial<int>(int(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  RPolynomial<int> f1 = p1.abs();
  RPolynomial<int> f2 = p2.abs();
  int f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  int F = ring_or_field<int>::gcd(f1c,f2c);
  RPolynomial<int> q,r; int M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    RPolynomial<int>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return RPolynomial<int>(F)*f1.abs();
}


// SPECIALIZING to :
 /*CGAL_KERNEL_MEDIUM_INLINE*/ 
void RPolynomial<double>::euclidean_div(
  const RPolynomial<double>& f, const RPolynomial<double>& g,
  RPolynomial<double>& q, RPolynomial<double>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = RPolynomial<double>(double(0)); }
  else { qd = rd-gd+1; q = RPolynomial<double>(size_t(qd)); }
  while ( rd >= gd ) {
    double S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}


 /*CGAL_KERNEL_MEDIUM_INLINE*/   
void RPolynomial<double>::pseudo_div(
  const RPolynomial<double>& f, const RPolynomial<double>& g, 
  RPolynomial<double>& q, RPolynomial<double>& r, double& D)
{
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = RPolynomial<double>(0); r = f; D = 1; 
    CGAL_postcondition(RPolynomial<double>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = RPolynomial<double>( size_t(delta) );
  double G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  RPolynomial<double> res = RPolynomial<double>(D)*f;
  TRACEN("  pseudo_div start "<<res<<" "<<qd<<" "<<q.degree());
  while (qd >= 0) {
    double F = res[rd]; // highest order coeff of res
    double t = F/G;     // ensured to be integer by multiplication of D
    q.coeff(qd) = t;    // store q coeff
    res.minus_offsetmult(g,t,qd); 
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res;
  CGAL_postcondition(RPolynomial<double>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}


 /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<double> RPolynomial<double>::gcd(
  const RPolynomial<double>& p1, const RPolynomial<double>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return RPolynomial<double>(double(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  RPolynomial<double> f1 = p1.abs();
  RPolynomial<double> f2 = p2.abs();
  double f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  double F = ring_or_field<double>::gcd(f1c,f2c);
  RPolynomial<double> q,r; double M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    RPolynomial<double>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return RPolynomial<double>(F)*f1.abs();
}


// SPECIALIZE_FUNCTION ORIGINAL
template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
void RPolynomial<NT>::euclidean_div(
  const RPolynomial<NT>& f, const RPolynomial<NT>& g,
  RPolynomial<NT>& q, RPolynomial<NT>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = RPolynomial<NT>(NT(0)); }
  else { qd = rd-gd+1; q = RPolynomial<NT>(size_t(qd)); }
  while ( rd >= gd ) {
    NT S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/   
void RPolynomial<NT>::pseudo_div(
  const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
  RPolynomial<NT>& q, RPolynomial<NT>& r, NT& D)
{
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = RPolynomial<NT>(0); r = f; D = 1; 
    CGAL_postcondition(RPolynomial<NT>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = RPolynomial<NT>( size_t(delta) );
  NT G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  RPolynomial<NT> res = RPolynomial<NT>(D)*f;
  TRACEN("  pseudo_div start "<<res<<" "<<qd<<" "<<q.degree());
  while (qd >= 0) {
    NT F = res[rd]; // highest order coeff of res
    NT t = F/G;     // ensured to be integer by multiplication of D
    q.coeff(qd) = t;    // store q coeff
    res.minus_offsetmult(g,t,qd); 
    if (res.is_zero()) break;
    rd = res.degree();
    qd = rd - gd;
  }
  r = res;
  CGAL_postcondition(RPolynomial<NT>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> RPolynomial<NT>::gcd(
  const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return RPolynomial<NT>(NT(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  RPolynomial<NT> f1 = p1.abs();
  RPolynomial<NT> f2 = p2.abs();
  NT f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  NT F = ring_or_field<NT>::gcd(f1c,f2c);
  RPolynomial<NT> q,r; NT M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    RPolynomial<NT>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return RPolynomial<NT>(F)*f1.abs();
}


// SPECIALIZE_IMPLEMENTATION(NT,int double) END

CGAL_END_NAMESPACE

#endif  // CGAL_RPOLYNOMIAL_H


