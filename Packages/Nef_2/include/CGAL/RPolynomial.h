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

//#define CGAL_NO_LEDA_HANDLE
//#define CGAL_DEBUG_HANDLE_REP

#include <CGAL/basic.h>
#include <vector>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/IO/io.h>
#include <CGAL/IO/Istream_iterator.h> 
#include <CGAL/IO/Ostream_iterator.h> 

#undef _DEBUG
#define _DEBUG 3
#include <CGAL/Nef_2/debug.h>

#ifdef _MSC_VER
#define MSC_HACK_ARGDECL ,char dummy1, char dummy2
#define MSC_HACK_ARGINS  ,'x','x'
#else
#define MSC_HACK_ARGDECL
#define MSC_HACK_ARGINS 
#define rp_gcd gcd
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
};

template <>
struct ring_or_field<double> {
  typedef field_with_div kind;
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

template <class Forward_iterator>
typename std::iterator_traits<Forward_iterator>::value_type 
gcd_of_range(Forward_iterator its, Forward_iterator ite)
/*{\Mfunc calculates the greates common divisor of the
set of numbers $\{ |*its|, |*++its|, \ldots, |*it| \}$ of type |NT|,
where |++it == ite| and |NT| is the value type of |Forward_iterator|. 
\precond there exists a pairwise gcd operation |NT gcd(NT,NT)| and 
|its!=ite|.}*/
{ CGAL_assertion(its!=ite);
  std::iterator_traits<Forward_iterator>::value_type res = *its++;
  for(; its!=ite; ++its) res = (*its==0 ? res : gcd(res, *its));
  if (res==0) res = 1;
  return res;
}



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


#ifdef _MSC_VER
/* no koenig lookup thus we put the following two completing
   functions into global space */
CGAL_END_NAMESPACE
#endif

int gcd(const int& a, const int& b)
{ if (a ==0 )
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

double gcd (const double&, const double&)
{ return 1.0; }

#ifdef _MSC_VER
CGAL_BEGIN_NAMESPACE
#endif

template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
  rp_gcd(const RPolynomial<NT>&, const RPolynomial<NT>&);

#ifndef _MSC_VER 
template<class NT> /*CGAL_KERNEL_INLINE*/ CGAL::Sign 
  sign(const RPolynomial<NT>& p);
#endif // collides with global CGAL sign

template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ 
  void pseudo_div(const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
         RPolynomial<NT>& q, RPolynomial<NT>& r, NT& D);

template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ 
  void euclidean_div(const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
         RPolynomial<NT>& q, RPolynomial<NT>& r);


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
  typedef typename std::vector<NT>::size_type Size_type;
  std::vector<NT> coeff;

  RPolynomial_rep() : coeff() {}
  RPolynomial_rep(const NT& n) : coeff(1) { coeff[0]=n; }
  RPolynomial_rep(const NT& n, const NT& m) : coeff(2)
    { coeff[0]=n; coeff[1]=m; }
  RPolynomial_rep(const NT& a, const NT& b, const NT& c) : coeff(3)
    { coeff[0]=a; coeff[1]=b; coeff[2]=c; }
  RPolynomial_rep(Size_type s) : coeff(s,NT(0)) {}

  template <class Forward_iterator>
  RPolynomial_rep(Forward_iterator first, Forward_iterator last
                  MSC_HACK_ARGDECL) : coeff(first,last) {}

  void reduce() 
  { typename std::vector<NT>::iterator it;
    while (*(--(it = coeff.end()))==NT(0) && coeff.size()>1)
    { coeff.pop_back(); }
  }

  friend class RPolynomial<pNT>;
  friend class RPolynomial<int>;
  friend class RPolynomial<double>;
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS  
         (std::istream&, RPolynomial<NT>&);


};

// SPECIALIZE_CLASS(NT,int double) START
// CLASS TEMPLATE NT: 
/*{\Msubst 
rp_gcd#gcd
typename iterator_traits<Forward_iterator>::value_type#NT
CGAL_NULL_TMPL_ARGS#
MSC_HACK_ARGDECL#
CGAL_SCOPE#
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

  typedef typename std::vector<NT>::const_iterator const_iterator;
  /*{\Mtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  typedef Handle_for< RPolynomial_rep<NT> > Base;
  typedef typename RPolynomial_rep<pNT>::Size_type Size_type;

  protected:
  void reduce() { ptr->reduce(); }
  std::vector<NT>& coeffs() { return ptr->coeff; }
  const std::vector<NT>& coeffs() const { return ptr->coeff; }
  RPolynomial(Size_type s) : Base( RPolynomial_rep<NT>(s) ) {}
  // creates a polynomial of degree s-1

  static NT _R; // for visualization only

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

  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last
    MSC_HACK_ARGDECL)

  /*{\Mcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<NT>(first,last MSC_HACK_ARGINS)) { reduce(); }

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
  /*{\Moperations 3 3 }*/

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

  protected: // accessing coefficients internally:
  NT& a(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared()); 
    CGAL_assertion(i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Mop a random access iterator pointing to $a_0$.}*/

  const_iterator end() const { return ptr->coeff.end(); }
  /*{\Mop a random access iterator pointing beyond $a_d$.}*/

  NT eval_at(const NT& R) const
  /*{\Mop evaluates the polynomial at |R|.}*/
  { CGAL_assertion( degree()>=0 );
    NT res = ptr->coeff[0];
    NT x = _R;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=_R; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Mop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { std::vector<NT>::const_iterator it = (ptr->coeff.end()); --it;
    if (*it < NT(0)) return (CGAL::NEGATIVE);
    if (*it > NT(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Mop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==0; }

  RPolynomial<NT> abs() const
  /*{\Mop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  NT content() const
  /*{\Mop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |NT| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }


  static void set_R(const NT& R) { _R = R; }

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
  /*{\Mbinopfunc implements polynomial division of |p1| and
  |p2|. if |p1 = p2 * p3| then |p2| is returned. The result
  is undefined if |p3| does not exist in |NT[x]|. 
  The correct division algorithm is chosen according to a traits class
  |ring_or_field<NT>| provided by the user.  If |ring_or_field<NT>::kind
  == ring_with_gcd| then the division is done by \emph{pseudo division}
  based on a |gcd| operation of |NT|.  If |ring_or_field<NT>::kind ==
  field_with_div| then the division is done by \emph{euclidean division}
  based on the division operation of the field |NT|.

  \textbf{Note} that |NT=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Mtext \headerline{Non member functions}}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<NT> rp_gcd CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2);
  /*{\Mfunc returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |NT=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |NT| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE pseudo_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
     RPolynomial<NT>& q, RPolynomial<NT>& r, NT& D);
  /*{\Mfunc implements division with remainder on polynomials of 
  the ring |NT[x]|: $D*f = g*q + r$.  \precond |NT| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |NT|.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE euclidean_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
     RPolynomial<NT>& q, RPolynomial<NT>& r);
  /*{\Mfunc implements division with remainder on polynomials of 
  the ring |NT[x]|: $f = g*q + r$.  \precond |NT| is a field, i.e.,
  there exists a division operation on |NT|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<NT>& p) ;


  RPolynomial<NT>& operator += (const RPolynomial<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<NT>& operator -= (const RPolynomial<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) -= p1[i];
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
    a(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const NT& num)
  { copy_on_write();
    a(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const NT& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const NT& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (NT)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const int& 
    
  RPolynomial<NT>& operator += (const int& num)
  { copy_on_write();
    a(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const int& num)
  { copy_on_write();
    a(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (NT)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const double& 
    
  RPolynomial<NT>& operator += (const double& num)
  { copy_on_write();
    a(0) += (NT)num; return *this; }

  RPolynomial<NT>& operator -= (const double& num)
  { copy_on_write();
    a(0) -= (NT)num; return *this; }

  RPolynomial<NT>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (NT)num; 
    reduce(); return *this; }

  RPolynomial<NT>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (NT)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<NT>& p, const NT& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<NT> s(Size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.a(i) = b*p[i-k];
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
rp_gcd#gcd
 iterator_traits<Forward_iterator>::value_type#int
CGAL_NULL_TMPL_ARGS#
MSC_HACK_ARGDECL#
CGAL_SCOPE#
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

  typedef  std::vector<int>::const_iterator const_iterator;
  /*{\Xtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  typedef Handle_for< RPolynomial_rep<int> > Base;
  typedef  RPolynomial_rep<int>::Size_type Size_type;

  protected:
  void reduce() { ptr->reduce(); }
  std::vector<int>& coeffs() { return ptr->coeff; }
  const std::vector<int>& coeffs() const { return ptr->coeff; }
  RPolynomial(Size_type s) : Base( RPolynomial_rep<int>(s) ) {}
  // creates a polynomial of degree s-1

  static int _R; // for visualization only

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

  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last
    MSC_HACK_ARGDECL)

  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<int>(first,last MSC_HACK_ARGINS)) { reduce(); }

  // KILL double START
  RPolynomial(double n) : Base(RPolynomial_rep<int>(int(n))) { reduce(); }
  RPolynomial(double n1, double n2) 
    : Base(RPolynomial_rep<int>(int(n1),int(n2))) { reduce(); }
  // KILL double END

  RPolynomial(const RPolynomial<int>& p) : Base(p) {}
  /*{\Xoperations 3 3 }*/

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

  protected: // accessing coefficients internally:
  int& a(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared()); 
    CGAL_assertion(i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Xop a random access iterator pointing to $a_0$.}*/

  const_iterator end() const { return ptr->coeff.end(); }
  /*{\Xop a random access iterator pointing beyond $a_d$.}*/

  int eval_at(const int& R) const
  /*{\Xop evaluates the polynomial at |R|.}*/
  { CGAL_assertion( degree()>=0 );
    int res = ptr->coeff[0];
    int x = _R;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=_R; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Xop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { std::vector<int>::const_iterator it = (ptr->coeff.end()); --it;
    if (*it < int(0)) return (CGAL::NEGATIVE);
    if (*it > int(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Xop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==0; }

  RPolynomial<int> abs() const
  /*{\Xop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  int content() const
  /*{\Xop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |int| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }


  static void set_R(const int& R) { _R = R; }

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
  /*{\Xbinopfunc implements polynomial division of |p1| and
  |p2|. if |p1 = p2 * p3| then |p2| is returned. The result
  is undefined if |p3| does not exist in |int[x]|. 
  The correct division algorithm is chosen according to a traits class
  |ring_or_field<int>| provided by the user.  If |ring_or_field<int>::kind
  == ring_with_gcd| then the division is done by \emph{pseudo division}
  based on a |gcd| operation of |int|.  If |ring_or_field<int>::kind ==
  field_with_div| then the division is done by \emph{euclidean division}
  based on the division operation of the field |int|.

  \textbf{Note} that |int=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Xtext \headerline{Non member functions}}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<int> rp_gcd CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<int>& p1, const RPolynomial<int>& p2);
  /*{\Xfunc returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |int=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |int| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE pseudo_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<int>& f, const RPolynomial<int>& g, 
     RPolynomial<int>& q, RPolynomial<int>& r, int& D);
  /*{\Xfunc implements division with remainder on polynomials of 
  the ring |int[x]|: $D*f = g*q + r$.  \precond |int| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |int|.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE euclidean_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<int>& f, const RPolynomial<int>& g, 
     RPolynomial<int>& q, RPolynomial<int>& r);
  /*{\Xfunc implements division with remainder on polynomials of 
  the ring |int[x]|: $f = g*q + r$.  \precond |int| is a field, i.e.,
  there exists a division operation on |int|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<int>& p) ;


  RPolynomial<int>& operator += (const RPolynomial<int>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<int>& operator -= (const RPolynomial<int>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) -= p1[i];
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
    a(0) += (int)num; return *this; }

  RPolynomial<int>& operator -= (const int& num)
  { copy_on_write();
    a(0) -= (int)num; return *this; }

  RPolynomial<int>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (int)num; 
    reduce(); return *this; }

  RPolynomial<int>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (int)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const double& 
    
  RPolynomial<int>& operator += (const double& num)
  { copy_on_write();
    a(0) += (int)num; return *this; }

  RPolynomial<int>& operator -= (const double& num)
  { copy_on_write();
    a(0) -= (int)num; return *this; }

  RPolynomial<int>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (int)num; 
    reduce(); return *this; }

  RPolynomial<int>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (int)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<int>& p, const int& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<int> s(Size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.a(i) = b*p[i-k];
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
rp_gcd#gcd
 iterator_traits<Forward_iterator>::value_type#double
CGAL_NULL_TMPL_ARGS#
MSC_HACK_ARGDECL#
CGAL_SCOPE#
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

  typedef  std::vector<double>::const_iterator const_iterator;
  /*{\Xtypemember a random access iterator for read-only access to the
  coefficient vector.}*/

  typedef Handle_for< RPolynomial_rep<double> > Base;
  typedef  RPolynomial_rep<double>::Size_type Size_type;

  protected:
  void reduce() { ptr->reduce(); }
  std::vector<double>& coeffs() { return ptr->coeff; }
  const std::vector<double>& coeffs() const { return ptr->coeff; }
  RPolynomial(Size_type s) : Base( RPolynomial_rep<double>(s) ) {}
  // creates a polynomial of degree s-1

  static double _R; // for visualization only

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

  template <class Forward_iterator>
  RPolynomial(Forward_iterator first, Forward_iterator last
    MSC_HACK_ARGDECL)

  /*{\Xcreate introduces a variable |\Mvar| of type |\Mname| representing
  the polynomial whose coefficients are determined by the iterator range,
  i.e. let $(a_0 = |*first|, a_1 = |*++first|, \ldots a_d = |*it|)$, 
  where |++it == last| then |\Mvar| stores the polynomial $a_1 + a_2 x + 
  \ldots a_d x^d$.}*/
    : Base(RPolynomial_rep<double>(first,last MSC_HACK_ARGINS)) { reduce(); }

  // KILL int START
  RPolynomial(int n) : Base(RPolynomial_rep<double>(double(n))) { reduce(); }
  RPolynomial(int n1, int n2)
    : Base(RPolynomial_rep<double>(double(n1),double(n2))) { reduce(); }
  // KILL int END

  RPolynomial(const RPolynomial<double>& p) : Base(p) {}
  /*{\Xoperations 3 3 }*/

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

  protected: // accessing coefficients internally:
  double& a(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared()); 
    CGAL_assertion(i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  const_iterator begin() const { return ptr->coeff.begin(); }
  /*{\Xop a random access iterator pointing to $a_0$.}*/

  const_iterator end() const { return ptr->coeff.end(); }
  /*{\Xop a random access iterator pointing beyond $a_d$.}*/

  double eval_at(const double& R) const
  /*{\Xop evaluates the polynomial at |R|.}*/
  { CGAL_assertion( degree()>=0 );
    double res = ptr->coeff[0];
    double x = _R;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=_R; }
    return res; 
  }

  CGAL::Sign sign() const
  /*{\Xop returns the sign of the limit process for $x \rightarrow \infty$
  (the sign of the leading coefficient).}*/
  { std::vector<double>::const_iterator it = (ptr->coeff.end()); --it;
    if (*it < double(0)) return (CGAL::NEGATIVE);
    if (*it > double(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  /*{\Xop returns true iff |\Mvar| is the zero polynomial.}*/
  { return degree()==0 && ptr->coeff[0]==0; }

  RPolynomial<double> abs() const
  /*{\Xop returns |-\Mvar| if |\Mvar.sign()==NEGATIVE| and |\Mvar| 
  otherwise.}*/
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  double content() const
  /*{\Xop returns the content of |\Mvar| (the gcd of its coefficients).
  \precond Requires |double| to provide a |gdc| operation.}*/
  { CGAL_assertion( degree()>=0 );
    return gcd_of_range(ptr->coeff.begin(),ptr->coeff.end());
  }


  static void set_R(const double& R) { _R = R; }

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
  /*{\Xbinopfunc implements polynomial division of |p1| and
  |p2|. if |p1 = p2 * p3| then |p2| is returned. The result
  is undefined if |p3| does not exist in |double[x]|. 
  The correct division algorithm is chosen according to a traits class
  |ring_or_field<double>| provided by the user.  If |ring_or_field<double>::kind
  == ring_with_gcd| then the division is done by \emph{pseudo division}
  based on a |gcd| operation of |double|.  If |ring_or_field<double>::kind ==
  field_with_div| then the division is done by \emph{euclidean division}
  based on the division operation of the field |double|.

  \textbf{Note} that |double=int| quickly leads to overflow
  errors when using this operation.}*/

  /*{\Xtext \headerline{Non member functions}}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  RPolynomial<double> rp_gcd CGAL_NULL_TMPL_ARGS 
  (const RPolynomial<double>& p1, const RPolynomial<double>& p2);
  /*{\Xfunc returns the greatest common divisor of |p1| and |p2|.
  \textbf{Note} that |double=int| quickly leads to overflow errors when
  using this operation.  \precond Requires |double| to be a unique
  factorization domain, i.e. to provide a |gdc| operation.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE pseudo_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<double>& f, const RPolynomial<double>& g, 
     RPolynomial<double>& q, RPolynomial<double>& r, double& D);
  /*{\Xfunc implements division with remainder on polynomials of 
  the ring |double[x]|: $D*f = g*q + r$.  \precond |double| is a unique
  factorization domain, i.e., there exists a |gcd| operation and an
  integral division operation on |double|.}*/

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ 

  void CGAL_SCOPE euclidean_div CGAL_NULL_TMPL_ARGS 
    (const RPolynomial<double>& f, const RPolynomial<double>& g, 
     RPolynomial<double>& q, RPolynomial<double>& r);
  /*{\Xfunc implements division with remainder on polynomials of 
  the ring |double[x]|: $f = g*q + r$.  \precond |double| is a field, i.e.,
  there exists a division operation on |double|.  }*/

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial<double>& p) ;


  RPolynomial<double>& operator += (const RPolynomial<double>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial<double>& operator -= (const RPolynomial<double>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) a(i) -= p1[i];
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
    a(0) += (double)num; return *this; }

  RPolynomial<double>& operator -= (const double& num)
  { copy_on_write();
    a(0) -= (double)num; return *this; }

  RPolynomial<double>& operator *= (const double& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (double)num; 
    reduce(); return *this; }

  RPolynomial<double>& operator /= (const double& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (double)num; 
    reduce(); return *this; }
// SPECIALIZING_MEMBERS FOR const int& 
    
  RPolynomial<double>& operator += (const int& num)
  { copy_on_write();
    a(0) += (double)num; return *this; }

  RPolynomial<double>& operator -= (const int& num)
  { copy_on_write();
    a(0) -= (double)num; return *this; }

  RPolynomial<double>& operator *= (const int& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) a(i) *= (double)num; 
    reduce(); return *this; }

  RPolynomial<double>& operator /= (const int& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) a(i) /= (double)num; 
    reduce(); return *this; }

  // SPECIALIZE_MEMBERS(int double) END
  //------------------------------------------------------------------

  void minus_offsetmult(const RPolynomial<double>& p, const double& b, int k)
  { CGAL_assertion(!ptr->is_shared());
    RPolynomial<double> s(Size_type(p.degree()+k+1)); // zero entries
    for (int i=k; i <= s.degree(); ++i) s.a(i) = b*p[i-k];
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

template <class NT> NT RPolynomial<NT>::_R;
int    RPolynomial<int>::_R;
double RPolynomial<double>::_R;

template <class NT> /*CGAL_KERNEL_INLINE*/ double to_double 
  (const RPolynomial<NT>& p) 
  { return (CGAL::to_double(p.eval_at(RPolynomial<NT>::_R))); }

template <class NT>  /*CGAL_KERNEL_INLINE*/ bool is_valid 
  (const RPolynomial<NT>& p) 
  { return (CGAL::is_valid(p[0])); }


template <class NT> /*CGAL_KERNEL_INLINE*/ bool is_finite 
  (const RPolynomial<NT>& p) 
  { return (CGAL::is_finite(p[0])); }

template <class NT> /*CGAL_KERNEL_INLINE*/ CGAL::io_Operator 
  io_tag(const RPolynomial<NT>&) 
  { return CGAL::io_Operator(); }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator - (const RPolynomial<NT>& p)
{
  CGAL_assertion(p.degree()>=0);
  RPolynomial<NT> res(p.coeffs().begin(),p.coeffs().end() MSC_HACK_ARGINS);
  std::vector<NT>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator + (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{ 
  typedef typename RPolynomial<NT>::Size_type Size_type;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial<NT>  p( (Size_type)(max + 1));
  for (i = 0; i <= min; ++i ) p.a(i) = p1[i]+p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.a(i)=p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.a(i)=p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator - (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{ 
  typedef typename RPolynomial<NT>::Size_type Size_type;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial<NT>  p( (Size_type)(max+1) );
  for (i = 0; i <= min; ++i ) p.a(i)=p1[i]-p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.a(i)= -p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.a(i)=  p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> operator * (const RPolynomial<NT>& p1, 
                            const RPolynomial<NT>& p2)
{
  typedef typename RPolynomial<NT>::Size_type Size_type;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  RPolynomial<NT>  p( (Size_type)(p1.degree()+p2.degree()+1)); 
  // initialize with zeros
  for (int i=0; i <= p1.degree(); ++i)
    for (int j=0; j <= p2.degree(); ++j)
    { p.a(i+j) += (p1[i]*p2[j]); }
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
void euclidean_div(const CGAL::RPolynomial<NT>& f,
                   const CGAL::RPolynomial<NT>& g,
                   RPolynomial<NT>& q, RPolynomial<NT>& r)
{
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = RPolynomial<NT>(size_t(0)); }
  else { qd = rd-gd+1; q = RPolynomial<NT>(size_t(qd)); }
  while ( rd >= gd ) {
    NT S = r[rd] / g[gd];
    qd = rd-gd;
    q.a(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial<NT> divop (const RPolynomial<NT>& p1, 
                       const RPolynomial<NT>& p2,
                       field_with_div)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  RPolynomial<NT> q,r;
  euclidean_div(p1,p2,q,r);
  CGAL_postcondition( (p2*q+r==p1) );
  return q;
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/   
void pseudo_div(const RPolynomial<NT>& f, const RPolynomial<NT>& g, 
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
    q.a(qd) = t;    // store q coeff
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
RPolynomial<NT> divop (const RPolynomial<NT>& p1, const RPolynomial<NT>& p2,
                       ring_with_gcd)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  RPolynomial<NT> q,r; NT D; 
  pseudo_div(p1,p2,q,r,D); 
  CGAL_postcondition( (p2*q+r==p1*RPolynomial<NT>(D)) );
  return q/=D;
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial<NT>
rp_gcd(const RPolynomial<NT>& p1, const RPolynomial<NT>& p2)
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
  NT F = gcd(f1c,f2c);
  RPolynomial<NT> q,r; NT M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return RPolynomial<NT>(F)*f1.abs();
}




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
        if (p[i]!=0) { os << " + "; print_monomial(os,p[i],i); }
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
        std::vector<NT> coeffs(d+1);
        for(i=0; i<=d; ++i) is >> coeffs[i];
        p = RPolynomial<NT>(coeffs.begin(),coeffs.end() MSC_HACK_ARGINS);
      }
      break;
    case CGAL::IO::BINARY :
      CGAL::read(is, d);
      if (d < 0) p = RPolynomial<NT>();
      else {
        std::vector<NT> coeffs(d+1);
        for(i=0; i<=d; ++i) 
        { CGAL::read(is,c); coeffs[i]=c; }
        p = RPolynomial<NT>(coeffs.begin(),coeffs.end() MSC_HACK_ARGINS);
      }
      break;
    default:
      CGAL_assertion_msg(0,"\nStream must be in ascii or binary mode\n");
      break;
  }
  return is;
}



CGAL_END_NAMESPACE

#endif  // CGAL_RPOLYNOMIAL_H


