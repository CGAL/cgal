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
// file          : include/CGAL/RPolynomial_MSC.h
// package       : Nef_2 
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// maintainer    : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Polynomials in one variable for MSC shit
// ============================================================================

#ifndef CGAL_RPOLYNOMIAL_MSC_H
#define CGAL_RPOLYNOMIAL_MSC_H

#include <CGAL/basic.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/IO/io.h>
#undef _DEBUG
#define _DEBUG 3
#include <CGAL/Nef_2/debug.h>

#include <CGAL/Nef_2/vector_MSC.h>
#define CGAL_SIMPLE_NEF_INTERFACE
#define SNIHACK ,char,char
#define SNIINST ,'c','c'

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
  static RT gcd(const RT& a, const RT& b) 
  { return 1.0; }
};


CGAL_BEGIN_NAMESPACE

template <class NT> class RPolynomial_rep_MSC;
template <class NT> class RPolynomial_MSC;

/*{\Mtext \headerline{Range template}}*/

template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
  operator - (const RPolynomial_MSC<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  RPolynomial_MSC<NT>
  operator + (const RPolynomial_MSC<NT>&, const RPolynomial_MSC<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  RPolynomial_MSC<NT>
  operator - (const RPolynomial_MSC<NT>&, const RPolynomial_MSC<NT>&);
template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
  operator * (const RPolynomial_MSC<NT>&, const RPolynomial_MSC<NT>&);
template <class NT> inline RPolynomial_MSC<NT>
  operator / (const RPolynomial_MSC<NT>&, const RPolynomial_MSC<NT>&);

template <class NT> /*CGAL_KERNEL_INLINE*/ double 
  to_double(const RPolynomial_MSC<NT>& p);
template <class NT>  /*CGAL_KERNEL_INLINE*/ bool 
  is_valid(const RPolynomial_MSC<NT>& p);
template <class NT> /*CGAL_KERNEL_INLINE*/ bool 
  is_finite(const RPolynomial_MSC<NT>& p);

template<class NT>  
  std::ostream& operator << (std::ostream& os, const RPolynomial_MSC<NT>& p);
template <class NT>  
  std::istream& operator >> (std::istream& is, RPolynomial_MSC<NT>& p);

template <class pNT> class RPolynomial_rep_MSC : public Ref_counted
{ 
  typedef pNT NT;
  typedef CGAL::vector_MSC<NT> Vector;
  typedef typename Vector::size_type      Size_type;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;
  Vector coeff;

  RPolynomial_rep_MSC() : coeff() {}
  RPolynomial_rep_MSC(const NT& n) : coeff(1) { coeff[0]=n; }
  RPolynomial_rep_MSC(const NT& n, const NT& m) : coeff(2)
    { coeff[0]=n; coeff[1]=m; }
  RPolynomial_rep_MSC(const NT& a, const NT& b, const NT& c) : coeff(3)
    { coeff[0]=a; coeff[1]=b; coeff[2]=c; }
  RPolynomial_rep_MSC(Size_type s) : coeff(s,NT(0)) {}

  template <class Forward_iterator>
  RPolynomial_rep_MSC(Forward_iterator first, Forward_iterator last SNIHACK) 
    : coeff() 
  { while (first!=last) coeff.push_back(*first++); }

  void reduce() 
  { while ( coeff.size()>1 && coeff.back()==NT(0) ) coeff.pop_back(); }

  friend class RPolynomial_MSC<pNT>;
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS  
         (std::istream&, RPolynomial_MSC<NT>&);

};

template <class pNT> class RPolynomial_MSC : 
  public Handle_for< RPolynomial_rep_MSC<pNT> >
{
public:
  typedef pNT NT;
  typedef Handle_for< RPolynomial_rep_MSC<NT> > Base;
  typedef RPolynomial_rep_MSC<NT> Rep;
  typedef typename Rep::Vector    Vector;
  typedef typename Rep::Size_type Size_type;
  typedef typename Rep::iterator  iterator;
  typedef typename Rep::const_iterator const_iterator;

protected:
  struct init_by_degree { init_by_degree(){} };
  void reduce() { ptr->reduce(); }
  Vector& coeffs() { return ptr->coeff; }
  const Vector& coeffs() const { return ptr->coeff; }
  RPolynomial_MSC(init_by_degree, Size_type s) : Base(Rep(s)) {}
  // creates a polynomial of degree s-1

public:
  static NT RR;

  RPolynomial_MSC() : Base(Rep()) {}
  RPolynomial_MSC(const NT& a0) : Base(Rep(a0)) { reduce(); }
  RPolynomial_MSC(const NT& a0, const NT& a1) : Base(Rep(a0,a1)) { reduce(); }
  RPolynomial_MSC(const NT& a0, const NT& a1,const NT& a2)
    : Base(Rep(a0,a1,a2)) { reduce(); }

  #define RPOL(I)\
  RPolynomial_MSC(I first, I last):Base(Rep(first,last SNIINST)){ reduce(); }
  RPOL(const NT*)
  #undef RPOL

  RPolynomial_MSC(const RPolynomial_MSC<NT>& p) : Base(p) {}

  int degree() const 
  { return ptr->coeff.size()-1; } 

  const NT& operator[](unsigned int i) const 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }

  const NT& operator[](unsigned int i) 
  { CGAL_assertion( i<(ptr->coeff.size()) );
    return ptr->coeff[i]; }

  protected: // accessing coefficients internally:
  NT& coeff(unsigned int i) 
  { CGAL_assertion(!ptr->is_shared()); 
    CGAL_assertion(i<(ptr->coeff.size()));
    return ptr->coeff[i]; 
  }
  public:

  const_iterator begin() const { return ptr->coeff.begin(); }
  const_iterator end() const { return ptr->coeff.end(); }

  NT eval_at(const NT& r) const
  { CGAL_assertion( degree()>=0 );
    NT res = ptr->coeff[0];
    NT x = r;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr->coeff[i]*x; x*=r; }
    return res; 
  }

  CGAL::Sign sign() const
  { const_iterator it = (ptr->coeff.end()); --it;
    if (*it < NT(0)) return (CGAL::NEGATIVE);
    if (*it > NT(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  { return degree()==0 && ptr->coeff[0]==NT(0); }

  RPolynomial_MSC<NT> abs() const
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  NT content() const
  { CGAL_assertion( degree()>=0 );
    iterator its=ptr->coeff.begin(),ite=ptr->coeff.end();
    NT res = *its++;
    for(; its!=ite; ++its) res = 
      (*its==0 ? res : ring_or_field<NT>::gcd(res, *its));
    if (res==0) res = 1;
    return res;
  }

  static void set_R(const NT& R) { RR = R; }

  friend  /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
    operator - CGAL_NULL_TMPL_ARGS  (const RPolynomial_MSC<NT>&);   
                          
  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
    operator + CGAL_NULL_TMPL_ARGS (const RPolynomial_MSC<NT>&, 
                                    const RPolynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
    operator - CGAL_NULL_TMPL_ARGS (const RPolynomial_MSC<NT>&, 
                                    const RPolynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>
    operator * CGAL_NULL_TMPL_ARGS (const RPolynomial_MSC<NT>&, 
                                    const RPolynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ RPolynomial_MSC<NT>  
    operator / CGAL_NULL_TMPL_ARGS (const RPolynomial_MSC<NT>&, 
				    const RPolynomial_MSC<NT>&);

static RPolynomial_MSC<NT> gcd
    (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return RPolynomial_MSC<NT>(NT(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  RPolynomial_MSC<NT> f1 = p1.abs();
  RPolynomial_MSC<NT> f2 = p2.abs();
  NT f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  NT F = ring_or_field<NT>::gcd(f1c,f2c);
  RPolynomial_MSC<NT> q,r; NT M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    RPolynomial_MSC<NT>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return RPolynomial_MSC<NT>(F)*f1.abs();
}


static void pseudo_div
    (const RPolynomial_MSC<NT>& f, const RPolynomial_MSC<NT>& g, 
     RPolynomial_MSC<NT>& q, RPolynomial_MSC<NT>& r, NT& D)
{ init_by_degree IBD;
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = RPolynomial_MSC<NT>(0); r = f; D = 1; 
    CGAL_postcondition(RPolynomial_MSC<NT>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = RPolynomial_MSC<NT>( IBD, Size_type(delta) );
  NT G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  RPolynomial_MSC<NT> res = RPolynomial_MSC<NT>(D)*f;
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
  CGAL_postcondition(RPolynomial_MSC<NT>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}



static void euclidean_div 
    (const RPolynomial_MSC<NT>& f, const RPolynomial_MSC<NT>& g, 
     RPolynomial_MSC<NT>& q, RPolynomial_MSC<NT>& r)
{
  init_by_degree IBD;
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = RPolynomial_MSC<NT>(IBD,Size_type(0)); }
  else { qd = rd-gd+1; q = RPolynomial_MSC<NT>(IBD,Size_type(qd)); }
  while ( rd >= gd ) {
    NT S = r[rd] / g[gd];
    qd = rd-gd;
    q.coeff(qd) += S;
    r.minus_offsetmult(g,S,qd);
    rd = r.degree();
  }
  CGAL_postcondition( f==q*g+r );
}
 

  friend /*CGAL_KERNEL_INLINE*/ double to_double
  CGAL_NULL_TMPL_ARGS (const RPolynomial_MSC<NT>& p);


  RPolynomial_MSC<NT>& operator += (const RPolynomial_MSC<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) += p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  RPolynomial_MSC<NT>& operator -= (const RPolynomial_MSC<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) -= p1[i];
    while (i<=p1.degree()) ptr->coeff.push_back(-p1[i++]);
    reduce(); return (*this); }

  RPolynomial_MSC<NT>& operator *= (const RPolynomial_MSC<NT>& p1)
  { (*this)=(*this)*p1; return (*this); }

  RPolynomial_MSC<NT>& operator /= (const RPolynomial_MSC<NT>& p1)
  { (*this)=(*this)/p1; return (*this); }

  RPolynomial_MSC<NT>& operator += (const NT& num)
  { copy_on_write();
    coeff(0) += num; return *this; }

  RPolynomial_MSC<NT>& operator -= (const NT& num)
  { copy_on_write();
    coeff(0) -= num; return *this; }

  RPolynomial_MSC<NT>& operator *= (const NT& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= num; 
    reduce(); return *this; }

  RPolynomial_MSC<NT>& operator /= (const NT& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= num; 
    reduce(); return *this; }
   
  void minus_offsetmult(const RPolynomial_MSC<NT>& p, const NT& b, int k)
  { CGAL_assertion(!ptr->is_shared()); init_by_degree IBD;
    RPolynomial_MSC<NT> s( IBD, Size_type(p.degree()+k+1) ); 
    for (int i=k; i <= s.degree(); ++i) s.coeff(i) = b*p[i-k];
    operator-=(s);
  }

};


template <class NT> NT RPolynomial_MSC<NT>::RR;

template <class NT> /*CGAL_KERNEL_INLINE*/ double to_double 
  (const RPolynomial_MSC<NT>& p) 
  { return (CGAL::to_double(p.eval_at(RPolynomial_MSC<NT>::RR))); }

template <class NT>  /*CGAL_KERNEL_INLINE*/ bool is_valid 
  (const RPolynomial_MSC<NT>& p) 
  { return (CGAL::is_valid(p[0])); }


template <class NT> /*CGAL_KERNEL_INLINE*/ bool is_finite 
  (const RPolynomial_MSC<NT>& p) 
  { return CGAL::is_finite(p[0]); }

template <class NT> /*CGAL_KERNEL_INLINE*/ CGAL::io_Operator 
  io_tag(const RPolynomial_MSC<NT>&) 
  { CGAL::io_Operator OP; return OP; }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> operator - (const RPolynomial_MSC<NT>& p)
{
  CGAL_assertion(p.degree()>=0);
  RPolynomial_MSC<NT> res(p.coeffs().begin(),p.coeffs().end());
  typename RPolynomial_MSC<NT>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> operator + (const RPolynomial_MSC<NT>& p1, 
                            const RPolynomial_MSC<NT>& p2)
{ 
  typedef typename RPolynomial_MSC<NT>::Size_type Size_type;
  typename RPolynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial_MSC<NT>  p( IBD, Size_type(max + 1) );
  for (i = 0; i <= min; ++i ) p.coeff(i) = p1[i]+p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)=p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> operator - (const RPolynomial_MSC<NT>& p1, 
                            const RPolynomial_MSC<NT>& p2)
{ 
  typedef typename RPolynomial_MSC<NT>::Size_type Size_type;
  typename RPolynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  RPolynomial_MSC<NT>  p( IBD, Size_type(max+1) );
  for (i = 0; i <= min; ++i ) p.coeff(i)=p1[i]-p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)= -p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=  p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> operator * (const RPolynomial_MSC<NT>& p1, 
                            const RPolynomial_MSC<NT>& p2)
{
  typedef typename RPolynomial_MSC<NT>::Size_type Size_type;
  typename RPolynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  RPolynomial_MSC<NT>  p( IBD, Size_type(p1.degree()+p2.degree()+1) ); 
  for (int i=0; i <= p1.degree(); ++i)
    for (int j=0; j <= p2.degree(); ++j)
    { p.coeff(i+j) += (p1[i]*p2[j]); }
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> divop (const RPolynomial_MSC<NT>& p1, 
                           const RPolynomial_MSC<NT>& p2,
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
  return RPolynomial_MSC<NT>(); // never reached
}


template <class NT> inline
RPolynomial_MSC<NT> operator / (const RPolynomial_MSC<NT>& p1, 
                                const RPolynomial_MSC<NT>& p2)
{ typedef typename ring_or_field<NT>::kind KIND;
  return divop(p1,p2,KIND()); }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> divop (const RPolynomial_MSC<NT>& p1, 
			   const RPolynomial_MSC<NT>& p2,
			   field_with_div)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  RPolynomial_MSC<NT> q,r;
  RPolynomial_MSC<NT>::euclidean_div(p1,p2,q,r);
  CGAL_postcondition( (p2*q+r==p1) );
  return q;
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
RPolynomial_MSC<NT> divop (const RPolynomial_MSC<NT>& p1, 
			   const RPolynomial_MSC<NT>& p2,
			   ring_with_gcd)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return RPolynomial_MSC<NT>(NT(0));
  RPolynomial_MSC<NT> q,r; NT D; 
  RPolynomial_MSC<NT>::pseudo_div(p1,p2,q,r,D); 
  CGAL_postcondition( (p2*q+r==p1*RPolynomial_MSC<NT>(D)) );
  return q/=D;
}


template <class NT> 
inline RPolynomial_MSC<NT> 
gcd(const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
{ return RPolynomial_MSC<NT>::gcd(p1,p2); }

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator == 
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator != 
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <  
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::NEGATIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <= 
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >  
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >= 
  (const RPolynomial_MSC<NT>& p1, const RPolynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::NEGATIVE ); }    

template <class NT> 
void print_monomial(std::ostream& os, const NT& n, int i)
{
  if (i==0) os << n;
  if (i==1) os << n << "R";
  if (i>1)  os << n << "R^" << i;
}


// I/O 
template <class NT>
std::ostream& operator << (std::ostream& os, const RPolynomial_MSC<NT>& p)
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
std::istream& operator >> (std::istream& is, RPolynomial_MSC<NT>& p)
{ 
  int i,d;
  NT c;
  switch( is.iword(CGAL::IO::mode) )
  { 
    case CGAL::IO::ASCII : 
      is >> d;
      if (d < 0) p = RPolynomial_MSC<NT>();
      else {
        typename RPolynomial_MSC<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) is >> coeffs[i];
        p = RPolynomial_MSC<NT>(coeffs.begin(),coeffs.end());
      }
      break;
    case CGAL::IO::BINARY :
      CGAL::read(is, d);
      if (d < 0) p = RPolynomial_MSC<NT>();
      else {
        typename RPolynomial_MSC<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) 
        { CGAL::read(is,c); coeffs[i]=c; }
        p = RPolynomial_MSC<NT>(coeffs.begin(),coeffs.end());
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


