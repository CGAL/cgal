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
// file          : include/CGAL/Nef_2/Polynomial_MSC.h
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

#ifndef CGAL_POLYNOMIAL_MSC_H
#define CGAL_POLYNOMIAL_MSC_H

#include <CGAL/basic.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/IO/io.h>
#undef _DEBUG
#define _DEBUG 3
#include <CGAL/Nef_2/debug.h>

#include <CGAL/Nef_2/vector_MSC.h>
#define CGAL_SIMPLE_NEF_INTERFACE
#define SNIHACK ,char,char
#define SNIINST ,'c','c'

CGAL_BEGIN_NAMESPACE

template <class NT> class Polynomial_rep_MSC;
template <class NT> class Polynomial_MSC;

/*{\Mtext \headerline{Range template}}*/

template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
  operator - (const Polynomial_MSC<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  Polynomial_MSC<NT>
  operator + (const Polynomial_MSC<NT>&, const Polynomial_MSC<NT>&);
template <class NT>  /*CGAL_KERNEL_MEDIUM_INLINE*/  Polynomial_MSC<NT>
  operator - (const Polynomial_MSC<NT>&, const Polynomial_MSC<NT>&);
template <class NT>   /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
  operator * (const Polynomial_MSC<NT>&, const Polynomial_MSC<NT>&);
template <class NT> inline Polynomial_MSC<NT>
  operator / (const Polynomial_MSC<NT>&, const Polynomial_MSC<NT>&);

template <class NT> /*CGAL_KERNEL_INLINE*/ double 
  to_double(const Polynomial_MSC<NT>& p);
template <class NT>  /*CGAL_KERNEL_INLINE*/ bool 
  is_valid(const Polynomial_MSC<NT>& p);
template <class NT> /*CGAL_KERNEL_INLINE*/ bool 
  is_finite(const Polynomial_MSC<NT>& p);

template<class NT>  
  std::ostream& operator << (std::ostream& os, const Polynomial_MSC<NT>& p);
template <class NT>  
  std::istream& operator >> (std::istream& is, Polynomial_MSC<NT>& p);

template <class pNT> class Polynomial_rep_MSC 
{ 
  typedef pNT NT;
  typedef CGAL::vector_MSC<NT> Vector;
  typedef typename Vector::size_type      Size_type;
  typedef typename Vector::iterator       iterator;
  typedef typename Vector::const_iterator const_iterator;
  Vector coeff;

  Polynomial_rep_MSC() : coeff() {}
  Polynomial_rep_MSC(const NT& n) : coeff(1) { coeff[0]=n; }
  Polynomial_rep_MSC(const NT& n, const NT& m) : coeff(2)
    { coeff[0]=n; coeff[1]=m; }
  Polynomial_rep_MSC(const NT& a, const NT& b, const NT& c) : coeff(3)
    { coeff[0]=a; coeff[1]=b; coeff[2]=c; }
  Polynomial_rep_MSC(Size_type s) : coeff(s,NT(0)) {}

  template <class Forward_iterator>
  Polynomial_rep_MSC(Forward_iterator first, Forward_iterator last SNIHACK) 
    : coeff() 
  { while (first!=last) coeff.push_back(*first++); }

  void reduce() 
  { while ( coeff.size()>1 && coeff.back()==NT(0) ) coeff.pop_back(); }

  friend class Polynomial_MSC<pNT>;
  friend std::istream& operator >> CGAL_NULL_TMPL_ARGS  
         (std::istream&, Polynomial_MSC<NT>&);

};

template <class pNT> class Polynomial_MSC : 
  public Handle_for< Polynomial_rep_MSC<pNT> >
{
public:
  typedef pNT NT;
  typedef Polynomial_rep_MSC<NT> Rep;
  typedef Handle_for< Rep >      Base;
  typedef typename Rep::Vector    Vector;
  typedef typename Rep::Size_type Size_type;
  typedef typename Rep::iterator  iterator;
  typedef typename Rep::const_iterator const_iterator;

protected:
  struct init_by_degree { init_by_degree(){} };
  void reduce() { ptr()->reduce(); }
  Vector& coeffs() { return ptr()->coeff; }
  const Vector& coeffs() const { return ptr()->coeff; }
  Polynomial_MSC(init_by_degree, Size_type s) : Base(Rep(s)) {}
  // creates a polynomial of degree s-1

public:
  static NT RR;

  Polynomial_MSC() : Base(Rep()) {}
  Polynomial_MSC(const NT& a0) : Base(Rep(a0)) { reduce(); }
  Polynomial_MSC(const NT& a0, const NT& a1) : Base(Rep(a0,a1)) { reduce(); }
  Polynomial_MSC(const NT& a0, const NT& a1,const NT& a2)
    : Base(Rep(a0,a1,a2)) { reduce(); }

  #define RPOL(I)\
  Polynomial_MSC(I first, I last):Base(Rep(first,last SNIINST)){ reduce(); }
  RPOL(const NT*)
  #undef RPOL

  Polynomial_MSC(const Polynomial_MSC<NT>& p) : Base(p) {}

  int degree() const 
  { return ptr()->coeff.size()-1; } 

  const NT& operator[](unsigned int i) const 
  { CGAL_assertion( i<(ptr()->coeff.size()) );
    return ptr()->coeff[i]; }

  const NT& operator[](unsigned int i) 
  { CGAL_assertion( i<(ptr()->coeff.size()) );
    return ptr()->coeff[i]; }

  protected: // accessing coefficients internally:
  NT& coeff(unsigned int i) 
  { CGAL_assertion(!is_shared()); 
    CGAL_assertion(i<(ptr()->coeff.size()));
    return ptr()->coeff[i]; 
  }
  public:

  const_iterator begin() const { return ptr()->coeff.begin(); }
  const_iterator end() const { return ptr()->coeff.end(); }

  NT eval_at(const NT& r) const
  { CGAL_assertion( degree()>=0 );
    NT res = ptr()->coeff[0];
    NT x = r;
    for(int i=1; i<=degree(); ++i) 
    { res += ptr()->coeff[i]*x; x*=r; }
    return res; 
  }

  CGAL::Sign sign() const
  { const_iterator it = (ptr()->coeff.end()); --it;
    if (*it < NT(0)) return (CGAL::NEGATIVE);
    if (*it > NT(0)) return (CGAL::POSITIVE);
    return CGAL::ZERO;
  }

  bool is_zero() const
  { return degree()==0 && ptr()->coeff[0]==NT(0); }

  Polynomial_MSC<NT> abs() const
  { if ( sign()==CGAL::NEGATIVE ) return -*this; return *this; }

  NT content() const
  { CGAL_assertion( degree()>=0 );
    const_iterator its=ptr()->coeff.begin(),ite=ptr()->coeff.end();
    NT res = *its++;
    for(; its!=ite; ++its) res = 
      (*its==0 ? res : CGAL_NTS gcd(res, *its));
    if (res==0) res = 1;
    return res;
  }

  static void set_R(const NT& R) { RR = R; }

  friend  /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
    operator - CGAL_NULL_TMPL_ARGS  (const Polynomial_MSC<NT>&);   
                          
  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
    operator + CGAL_NULL_TMPL_ARGS (const Polynomial_MSC<NT>&, 
                                    const Polynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
    operator - CGAL_NULL_TMPL_ARGS (const Polynomial_MSC<NT>&, 
                                    const Polynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>
    operator * CGAL_NULL_TMPL_ARGS (const Polynomial_MSC<NT>&, 
                                    const Polynomial_MSC<NT>&);

  friend /*CGAL_KERNEL_MEDIUM_INLINE*/ Polynomial_MSC<NT>  
    operator / CGAL_NULL_TMPL_ARGS (const Polynomial_MSC<NT>&, 
				    const Polynomial_MSC<NT>&);

static Polynomial_MSC<NT> gcd
    (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
{ TRACEN("gcd("<<p1<<" , "<<p2<<")");
  if ( p1.is_zero() )
    if ( p2.is_zero() ) return Polynomial_MSC<NT>(NT(1));
    else return p2.abs();
  if ( p2.is_zero() )
    return p1.abs();

  Polynomial_MSC<NT> f1 = p1.abs();
  Polynomial_MSC<NT> f2 = p2.abs();
  NT f1c = f1.content(), f2c = f2.content();
  f1 /= f1c; f2 /= f2c;
  NT F = CGAL_NTS gcd(f1c,f2c);
  Polynomial_MSC<NT> q,r; NT M=1,D;
  bool first = true;
  while ( ! f2.is_zero() ) { 
    Polynomial_MSC<NT>::pseudo_div(f1,f2,q,r,D);
    if (!first) M*=D;
    TRACEV(f1);TRACEV(f2);TRACEV(q);TRACEV(r);TRACEV(M);
    r /= r.content();
    f1=f2; f2=r;
    first=false;
  }
  TRACEV(f1.content());
  return Polynomial_MSC<NT>(F)*f1.abs();
}


static void pseudo_div
    (const Polynomial_MSC<NT>& f, const Polynomial_MSC<NT>& g, 
     Polynomial_MSC<NT>& q, Polynomial_MSC<NT>& r, NT& D)
{ init_by_degree IBD;
  TRACEN("pseudo_div "<<f<<" , "<< g);
  int fd=f.degree(), gd=g.degree();
  if ( fd<gd ) 
  { q = Polynomial_MSC<NT>(0); r = f; D = 1; 
    CGAL_postcondition(Polynomial_MSC<NT>(D)*f==q*g+r); return; 
  }
  // now we know fd >= gd and f>=g
  int qd=fd-gd, delta=qd+1, rd=fd;
  q = Polynomial_MSC<NT>( IBD, Size_type(delta) );
  NT G = g[gd]; // highest order coeff of g
  D = G; while (--delta) D*=G; // D = G^delta
  Polynomial_MSC<NT> res = Polynomial_MSC<NT>(D)*f;
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
  CGAL_postcondition(Polynomial_MSC<NT>(D)*f==q*g+r);
  TRACEN("  returning "<<q<<", "<<r<<", "<< D);
}



static void euclidean_div 
    (const Polynomial_MSC<NT>& f, const Polynomial_MSC<NT>& g, 
     Polynomial_MSC<NT>& q, Polynomial_MSC<NT>& r)
{
  init_by_degree IBD;
  r = f; r.copy_on_write();
  int rd=r.degree(), gd=g.degree(), qd(0);
  if ( rd < gd ) { q = Polynomial_MSC<NT>(IBD,Size_type(0)); }
  else { qd = rd-gd+1; q = Polynomial_MSC<NT>(IBD,Size_type(qd)); }
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
  CGAL_NULL_TMPL_ARGS (const Polynomial_MSC<NT>& p);


  Polynomial_MSC<NT>& operator += (const Polynomial_MSC<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) += p1[i];
    while (i<=p1.degree()) ptr()->coeff.push_back(p1[i++]);
    reduce(); return (*this); }

  Polynomial_MSC<NT>& operator -= (const Polynomial_MSC<NT>& p1)
  { copy_on_write();
    int d = std::min(degree(),p1.degree()), i;
    for(i=0; i<=d; ++i) coeff(i) -= p1[i];
    while (i<=p1.degree()) ptr()->coeff.push_back(-p1[i++]);
    reduce(); return (*this); }

  Polynomial_MSC<NT>& operator *= (const Polynomial_MSC<NT>& p1)
  { (*this)=(*this)*p1; return (*this); }

  Polynomial_MSC<NT>& operator /= (const Polynomial_MSC<NT>& p1)
  { (*this)=(*this)/p1; return (*this); }

  Polynomial_MSC<NT>& operator += (const NT& num)
  { copy_on_write();
    coeff(0) += num; return *this; }

  Polynomial_MSC<NT>& operator -= (const NT& num)
  { copy_on_write();
    coeff(0) -= num; return *this; }

  Polynomial_MSC<NT>& operator *= (const NT& num)
  { copy_on_write();
    for(int i=0; i<=degree(); ++i) coeff(i) *= num; 
    reduce(); return *this; }

  Polynomial_MSC<NT>& operator /= (const NT& num)
  { copy_on_write(); CGAL_assertion(num!=0);
    for(int i=0; i<=degree(); ++i) coeff(i) /= num; 
    reduce(); return *this; }
   
  void minus_offsetmult(const Polynomial_MSC<NT>& p, const NT& b, int k)
  { CGAL_assertion(!is_shared()); init_by_degree IBD;
    Polynomial_MSC<NT> s( IBD, Size_type(p.degree()+k+1) ); 
    for (int i=k; i <= s.degree(); ++i) s.coeff(i) = b*p[i-k];
    operator-=(s);
  }

};


template <class NT> NT Polynomial_MSC<NT>::RR;

template <class NT> /*CGAL_KERNEL_INLINE*/ double to_double 
  (const Polynomial_MSC<NT>& p) 
  { return (CGAL::to_double(p.eval_at(Polynomial_MSC<NT>::RR))); }

template <class NT>  /*CGAL_KERNEL_INLINE*/ bool is_valid 
  (const Polynomial_MSC<NT>& p) 
  { return (CGAL::is_valid(p[0])); }


template <class NT> /*CGAL_KERNEL_INLINE*/ bool is_finite 
  (const Polynomial_MSC<NT>& p) 
  { return CGAL::is_finite(p[0]); }

template <class NT> /*CGAL_KERNEL_INLINE*/ CGAL::io_Operator 
  io_tag(const Polynomial_MSC<NT>&) 
  { CGAL::io_Operator OP; return OP; }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> operator - (const Polynomial_MSC<NT>& p)
{
  CGAL_assertion(p.degree()>=0);
  Polynomial_MSC<NT> res(p.coeffs().begin(),p.coeffs().end());
  typename Polynomial_MSC<NT>::iterator it, ite=res.coeffs().end();
  for(it=res.coeffs().begin(); it!=ite; ++it) *it = -*it;
  return res;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> operator + (const Polynomial_MSC<NT>& p1, 
                            const Polynomial_MSC<NT>& p2)
{ 
  typedef typename Polynomial_MSC<NT>::Size_type Size_type;
  typename Polynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  Polynomial_MSC<NT>  p( IBD, Size_type(max + 1) );
  for (i = 0; i <= min; ++i ) p.coeff(i) = p1[i]+p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)=p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> operator - (const Polynomial_MSC<NT>& p1, 
                            const Polynomial_MSC<NT>& p2)
{ 
  typedef typename Polynomial_MSC<NT>::Size_type Size_type;
  typename Polynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  bool p1d_smaller_p2d = p1.degree() < p2.degree();
  int min,max,i;
  if (p1d_smaller_p2d) { min = p1.degree(); max = p2.degree(); }
  else                 { max = p1.degree(); min = p2.degree(); }
  Polynomial_MSC<NT>  p( IBD, Size_type(max+1) );
  for (i = 0; i <= min; ++i ) p.coeff(i)=p1[i]-p2[i];
  if (p1d_smaller_p2d)  for (; i <= max; ++i ) p.coeff(i)= -p2[i];
  else /* p1d >= p2d */ for (; i <= max; ++i ) p.coeff(i)=  p1[i];
  p.reduce();
  return p;
}

template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> operator * (const Polynomial_MSC<NT>& p1, 
                            const Polynomial_MSC<NT>& p2)
{
  typedef typename Polynomial_MSC<NT>::Size_type Size_type;
  typename Polynomial_MSC<NT>::init_by_degree IBD;
  CGAL_assertion(p1.degree()>=0);
  CGAL_assertion(p2.degree()>=0);
  Polynomial_MSC<NT>  p( IBD, Size_type(p1.degree()+p2.degree()+1) ); 
  for (int i=0; i <= p1.degree(); ++i)
    for (int j=0; j <= p2.degree(); ++j)
    { p.coeff(i+j) += (p1[i]*p2[j]); }
  p.reduce();
  return p;
}


template <class NT> inline
Polynomial_MSC<NT> operator / (const Polynomial_MSC<NT>& p1, 
                                const Polynomial_MSC<NT>& p2)
{ typedef typename Number_type_traits<NT>::Has_gcd HAS_GCD;
  return divop(p1,p2,HAS_GCD()); }


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> divop (const Polynomial_MSC<NT>& p1, 
			   const Polynomial_MSC<NT>& p2,
			   Tag_false)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return 0;
  Polynomial_MSC<NT> q,r;
  Polynomial_MSC<NT>::euclidean_div(p1,p2,q,r);
  CGAL_postcondition( (p2*q+r==p1) );
  return q;
}


template <class NT> /*CGAL_KERNEL_MEDIUM_INLINE*/ 
Polynomial_MSC<NT> divop (const Polynomial_MSC<NT>& p1, 
			   const Polynomial_MSC<NT>& p2,
			   Tag_true)
{ CGAL_assertion(!p2.is_zero());
  if (p1.is_zero()) return Polynomial_MSC<NT>(NT(0));
  Polynomial_MSC<NT> q,r; NT D; 
  Polynomial_MSC<NT>::pseudo_div(p1,p2,q,r,D); 
  CGAL_postcondition( (p2*q+r==p1*Polynomial_MSC<NT>(D)) );
  return q/=D;
}


template <class NT> 
inline Polynomial_MSC<NT> 
gcd(const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
{ return Polynomial_MSC<NT>::gcd(p1,p2); }

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator == 
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator != 
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::ZERO ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <  
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::NEGATIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator <= 
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() != CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >  
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
  { return ( (p1-p2).sign() == CGAL::POSITIVE ); }    

template <class NT> /*CGAL_KERNEL_INLINE*/ bool operator >= 
  (const Polynomial_MSC<NT>& p1, const Polynomial_MSC<NT>& p2)
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
std::ostream& operator << (std::ostream& os, const Polynomial_MSC<NT>& p)
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
      os << "Polynomial(" << p.degree() << ", ";
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
std::istream& operator >> (std::istream& is, Polynomial_MSC<NT>& p)
{ 
  int i,d;
  NT c;
  switch( is.iword(CGAL::IO::mode) )
  { 
    case CGAL::IO::ASCII : 
      is >> d;
      if (d < 0) p = Polynomial_MSC<NT>();
      else {
        typename Polynomial_MSC<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) is >> coeffs[i];
        p = Polynomial_MSC<NT>(coeffs.begin(),coeffs.end());
      }
      break;
    case CGAL::IO::BINARY :
      CGAL::read(is, d);
      if (d < 0) p = Polynomial_MSC<NT>();
      else {
        typename Polynomial_MSC<NT>::Vector coeffs(d+1);
        for(i=0; i<=d; ++i) 
        { CGAL::read(is,c); coeffs[i]=c; }
        p = Polynomial_MSC<NT>(coeffs.begin(),coeffs.end());
      }
      break;
    default:
      CGAL_assertion_msg(0,"\nStream must be in ascii or binary mode\n");
      break;
  }
  return is;
}


CGAL_END_NAMESPACE

#endif  // CGAL_POLYNOMIAL_H


